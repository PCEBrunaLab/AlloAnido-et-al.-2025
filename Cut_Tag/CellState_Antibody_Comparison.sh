#!/bin/bash
#SBATCH --job-name=CellState_Antibody_Comparison
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%j.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%j.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL

# Load modules
module use /opt/software/easybuild/modules/all/
module load Mamba
source ~/.bashrc
conda activate deeptools_fixed
module load SAMtools/1.11
module load bedtools/2.29.2

# =============================
# Parameters
# =============================
projPath="/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise"
OUT_DIR="${projPath}/heatmap/cellstate_condition_antibody"
GENES_BED="${projPath}/genes.bed"

# Cell state gene files
ADRN_GENE_FILE="${projPath}/Reference/vanGroningen_2017_ADRN.csv"
MES_GENE_FILE="${projPath}/Reference/vanGroningen_2017_MES.csv"
INT_GENE_FILE="${projPath}/Reference/Int.csv"

# Distance for Homer merged peaks
distance=7500

# Maximum distance for peak-to-gene assignment
MAX_DISTANCE=5000

# Markers and conditions
MARKERS=("H3K4me3" "KDM5A" "KDM5B")
CONDITIONS=("Untreated" "D7_Cis" "D14_Cis" "D14_C70")

# =============================
# Setup
# =============================
mkdir -p "${OUT_DIR}"

echo "================================================================"
echo "Cell State Condition-Antibody Comparison Analysis"
echo "Comparing H3K4me3, KDM5A, KDM5B at same positions per condition"
echo "Cell states: ADRN, MES, Int"
echo "Peak-to-gene assignment window: ${MAX_DISTANCE} bp"
echo "================================================================"

# Validate inputs
if [ ! -f "${GENES_BED}" ]; then
    echo "ERROR: genes.bed file not found at ${GENES_BED}"
    exit 1
fi

# =============================
# Process cell state gene files
# =============================
echo "Processing cell state gene lists..."

declare -A GENE_LISTS

# Process ADRN genes
if [ -f "${ADRN_GENE_FILE}" ]; then
    ADRN_LIST="${OUT_DIR}/ADRN_genes.txt"
    if head -1 "${ADRN_GENE_FILE}" | grep -q "Gene"; then
        tail -n +2 "${ADRN_GENE_FILE}" | awk -F',' '{print $1}' | tr -d ' \r' | sort -u > "${ADRN_LIST}"
    else
        awk -F',' '{print $1}' "${ADRN_GENE_FILE}" | tr -d ' \r' | sort -u > "${ADRN_LIST}"
    fi
    ADRN_COUNT=$(wc -l < "${ADRN_LIST}")
    echo "  ADRN: ${ADRN_COUNT} genes"
    GENE_LISTS["ADRN"]="${ADRN_LIST}"
fi

# Process MES genes
if [ -f "${MES_GENE_FILE}" ]; then
    MES_LIST="${OUT_DIR}/MES_genes.txt"
    if head -1 "${MES_GENE_FILE}" | grep -q "Gene"; then
        tail -n +2 "${MES_GENE_FILE}" | awk -F',' '{print $1}' | tr -d ' \r' | sort -u > "${MES_LIST}"
    else
        awk -F',' '{print $1}' "${MES_GENE_FILE}" | tr -d ' \r' | sort -u > "${MES_LIST}"
    fi
    MES_COUNT=$(wc -l < "${MES_LIST}")
    echo "  MES: ${MES_COUNT} genes"
    GENE_LISTS["MES"]="${MES_LIST}"
fi

# Process Int genes
if [ -f "${INT_GENE_FILE}" ]; then
    INT_LIST="${OUT_DIR}/Int_genes.txt"
    if head -1 "${INT_GENE_FILE}" | grep -q "Gene"; then
        tail -n +2 "${INT_GENE_FILE}" | awk -F',' '{print $1}' | tr -d ' \r' | sort -u > "${INT_LIST}"
    else
        awk -F',' '{print $1}' "${INT_GENE_FILE}" | tr -d ' \r' | sort -u > "${INT_LIST}"
    fi
    INT_COUNT=$(wc -l < "${INT_LIST}")
    echo "  Int: ${INT_COUNT} genes"
    GENE_LISTS["Int"]="${INT_LIST}"
fi

# =============================
# Create global union peak file for all antibodies
# =============================
UNION_PEAK_FILE="${OUT_DIR}/union_peaks_all_antibodies.bed"
echo ""
echo "Creating union peak file from all antibodies..."

> ${UNION_PEAK_FILE}.tmp
for hist in "${MARKERS[@]}"; do
    HOMER_PEAK="${projPath}/Homer/mergePeaks/human/Untreated_control/${hist}/${distance}/mergedPeakFile.bed"
    if [ -f "$HOMER_PEAK" ]; then
        cat "$HOMER_PEAK" >> ${UNION_PEAK_FILE}.tmp
    else
        echo "  Warning: Homer peaks not found for ${hist}, trying MACS2 peaks..."
        cat ${projPath}/peakCalling/MACS2/human/Untreated_control/*/${hist}/rep*/macs2_peak_q0.1_peaks.narrowPeak.bed 2>/dev/null >> ${UNION_PEAK_FILE}.tmp
    fi
done

if [ -s "${UNION_PEAK_FILE}.tmp" ]; then
    sort -k1,1 -k2,2n ${UNION_PEAK_FILE}.tmp | \
        bedtools merge -i - | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"union_peak_"NR,0,"."}' > ${UNION_PEAK_FILE}
    rm ${UNION_PEAK_FILE}.tmp
    GLOBAL_PEAK_COUNT=$(wc -l < ${UNION_PEAK_FILE})
    echo "Global union peak file created: ${GLOBAL_PEAK_COUNT} regions"
else
    echo "ERROR: No peaks found for any antibody"
    exit 1
fi

# =============================
# Process each cell state
# =============================
for state in "${!GENE_LISTS[@]}"; do
    echo ""
    echo "================================================================"
    echo "Processing cell state: ${state}"
    echo "================================================================"
    
    GENE_LIST="${GENE_LISTS[${state}]}"
    STATE_DIR="${OUT_DIR}/${state}"
    mkdir -p "${STATE_DIR}"
    
    # Create gene subset BED
    SUBSET_BED="${STATE_DIR}/genes_subset.bed"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]; next} ($4 in a)' \
        "${GENE_LIST}" "${GENES_BED}" | \
        sort -k1,1 -k2,2n | uniq > "${SUBSET_BED}"
    
    NUM_GENES=$(wc -l < "${SUBSET_BED}" 2>/dev/null || echo "0")
    if [ "${NUM_GENES}" -eq 0 ]; then
        echo "WARNING: No matching genes found for ${state}. Skipping."
        continue
    fi
    
    echo "Matched genes: ${NUM_GENES}"
    
    # =============================
    # Create TSS BED file
    # =============================
    TSS_BED="${STATE_DIR}/tss.bed"
    awk 'BEGIN{OFS="\t"} {
        if(NF<6) next;
        if($6=="+") tss=$2; 
        else if($6=="-") tss=$3-1; 
        else next;
        if(tss<0) tss=0;
        print $1, tss, tss+1, $4"_TSS", 0, $6
    }' "${SUBSET_BED}" | sort -k1,1 -k2,2n | uniq > "${TSS_BED}"
    
    echo "Created $(wc -l < "${TSS_BED}") TSS regions"
    
    # =============================
    # Create cell state-specific peak subset from union peaks
    # =============================
    echo "Filtering union peaks for ${state} genes..."
    
    # Check chromosome format consistency
    GENE_CHR=$(head -1 "${SUBSET_BED}" 2>/dev/null | cut -f1 | grep -c "^chr" || echo "0")
    PEAK_CHR=$(head -1 "${UNION_PEAK_FILE}" 2>/dev/null | cut -f1 | grep -c "^chr" || echo "0")
    
    PEAK_CLEAN="${STATE_DIR}/union_peaks_clean.bed"
    
    if [ "${GENE_CHR}" != "${PEAK_CHR}" ]; then
        echo "  Fixing chromosome naming mismatch..."
        if [ "${GENE_CHR}" -eq 1 ]; then
            # Genes have chr, add to peaks if missing
            awk 'BEGIN{OFS="\t"} {
                if(NF<3) next;
                if($2<0) $2=0;
                if($3<=$2) $3=$2+1;
                if($1 !~ /^chr/) $1="chr"$1;
                print $1,$2,$3,(NF>=4 ? $4 : "peak_"NR),0,"."
            }' "${UNION_PEAK_FILE}" | sort -k1,1 -k2,2n > "${PEAK_CLEAN}"
        else
            # Genes don't have chr, remove from peaks if present
            awk 'BEGIN{OFS="\t"} {
                if(NF<3) next;
                if($2<0) $2=0;
                if($3<=$2) $3=$2+1;
                gsub(/^chr/, "", $1);
                print $1,$2,$3,(NF>=4 ? $4 : "peak_"NR),0,"."
            }' "${UNION_PEAK_FILE}" | sort -k1,1 -k2,2n > "${PEAK_CLEAN}"
        fi
    else
        # No mismatch, just clean the file
        awk 'BEGIN{OFS="\t"} {
            if(NF<3) next;
            if($2<0) $2=0;
            if($3<=$2) $3=$2+1;
            print $1,$2,$3,(NF>=4 ? $4 : "peak_"NR),0,"."
        }' "${UNION_PEAK_FILE}" | sort -k1,1 -k2,2n > "${PEAK_CLEAN}"
    fi
    
    # Filter peaks near genes
    FILTERED_PEAKS="${STATE_DIR}/gene_associated_union_peaks.bed"
    
    bedtools window -w ${MAX_DISTANCE} -a "${PEAK_CLEAN}" -b "${SUBSET_BED}" | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_near_"$10,0,"."}' | \
        cut -f1-3 | sort -u | \
        awk 'BEGIN{OFS="\t"} {print $1,$2,$3,"filtered_peak_"NR,0,"."}' > "${FILTERED_PEAKS}"
    
    if [ ! -s "${FILTERED_PEAKS}" ]; then
        echo "  No peaks found within ${MAX_DISTANCE} bp of ${state} genes"
        FILTERED_PEAKS=""
    else
        PEAK_COUNT=$(wc -l < "${FILTERED_PEAKS}")
        echo "  Found ${PEAK_COUNT} union peaks associated with ${state} genes"
    fi
    
    # =============================
    # Process each condition
    # =============================
    for condition in "${CONDITIONS[@]}"; do
        echo ""
        echo "  Processing ${condition} for ${state} cell state..."
        
        # Build BigWig file list for all antibodies at this condition
        BW_FILES=""
        LABELS=""
        
        for hist in "${MARKERS[@]}"; do
            bw_file="${projPath}/alignment/bigwig/human/${condition}/${hist}/normalized.bw"
            if [ -f "$bw_file" ]; then
                BW_FILES="${BW_FILES} ${bw_file}"
                LABELS="${LABELS} ${hist}"
            else
                echo "    Warning: Missing BigWig for ${condition} ${hist}"
            fi
        done
        
        if [ -z "${BW_FILES}" ]; then
            echo "    ERROR: No BigWig files found for ${condition}"
            continue
        fi
        
        # Convert labels to array
        read -ra LABEL_ARRAY <<< "${LABELS}"
        
        # =============================
        # TSS-centered analysis
        # =============================
        echo "    Creating TSS-centered heatmap..."
        TSS_MATRIX="${STATE_DIR}/${condition}_TSS_matrix.gz"
        
        computeMatrix reference-point \
            -S ${BW_FILES} \
            -R "${TSS_BED}" \
            --referencePoint TSS \
            -a 3000 -b 3000 \
            --skipZeros \
            --missingDataAsZero \
            -o "${TSS_MATRIX}" \
            --samplesLabel ${LABEL_ARRAY[@]} \
            -p 8
        
        if [ -f "${TSS_MATRIX}" ]; then
            # TSS heatmap
            plotHeatmap \
                -m "${TSS_MATRIX}" \
                -out "${STATE_DIR}/${condition}_${state}_TSS_heatmap.pdf" \
                --colorMap RdBu_r \
                --whatToShow 'plot, heatmap and colorbar' \
                --heatmapHeight 12 \
                --heatmapWidth 8 \
                --plotTitle "${state} Cell State - ${condition} TSS Signal" \
                --xAxisLabel "Distance from TSS (bp)" \
                --yAxisLabel "Genes (n = ${NUM_GENES})" \
                --refPointLabel "TSS" \
                --sortRegions descend \
                --sortUsing mean \
                --zMin 0 --zMax auto
            
            # TSS profile plot
            plotProfile \
                -m "${TSS_MATRIX}" \
                -out "${STATE_DIR}/${condition}_${state}_TSS_profile.pdf" \
                --plotTitle "${state} Cell State - ${condition} Average TSS Signal" \
                --yAxisLabel "Normalized read density" \
                --xAxisLabel "Distance from TSS (bp)" \
                --refPointLabel "TSS" \
                --colors "#E41A1C" "#4DAF4A" "#984EA3" \
                --plotHeight 6 \
                --plotWidth 8 \
                --legendLocation upper-right \
                --perGroup
            
            echo "      ✓ TSS heatmap and profile created"
        fi
        
        # =============================
        # Gene body-scaled analysis
        # =============================
        echo "    Creating gene body heatmap..."
        GENE_BODY_MATRIX="${STATE_DIR}/${condition}_gene_body_matrix.gz"
        
        computeMatrix scale-regions \
            -S ${BW_FILES} \
            -R "${SUBSET_BED}" \
            --beforeRegionStartLength 3000 \
            --regionBodyLength 5000 \
            --afterRegionStartLength 3000 \
            --skipZeros \
            --missingDataAsZero \
            -o "${GENE_BODY_MATRIX}" \
            --samplesLabel ${LABEL_ARRAY[@]} \
            -p 8
        
        if [ -f "${GENE_BODY_MATRIX}" ]; then
            # Gene body heatmap
            plotHeatmap \
                -m "${GENE_BODY_MATRIX}" \
                -out "${STATE_DIR}/${condition}_${state}_gene_body_heatmap.pdf" \
                --colorMap RdBu_r \
                --whatToShow 'plot, heatmap and colorbar' \
                --heatmapHeight 12 \
                --heatmapWidth 8 \
                --plotTitle "${state} Cell State - ${condition} Gene Body Signal" \
                --xAxisLabel "Gene Body (scaled to 5kb)" \
                --yAxisLabel "Genes (n = ${NUM_GENES})" \
                --startLabel "Gene Start" \
                --endLabel "Gene End" \
                --sortRegions descend \
                --sortUsing mean \
                --zMin 0 --zMax auto
            
            # Gene body profile plot
            plotProfile \
                -m "${GENE_BODY_MATRIX}" \
                -out "${STATE_DIR}/${condition}_${state}_gene_body_profile.pdf" \
                --plotTitle "${state} Cell State - ${condition} Average Gene Body Signal" \
                --yAxisLabel "Normalized read density" \
                --startLabel "Gene Start" \
                --endLabel "Gene End" \
                --colors "#E41A1C" "#4DAF4A" "#984EA3" \
                --plotHeight 6 \
                --plotWidth 8 \
                --legendLocation upper-right \
                --perGroup
            
            echo "      ✓ Gene body heatmap and profile created"
        fi
        
        # =============================
        # Peak-centered analysis (if peaks exist)
        # =============================
        if [ -s "${FILTERED_PEAKS}" ]; then
            echo "    Creating peak-centered heatmap..."
            PEAK_MATRIX="${STATE_DIR}/${condition}_peaks_matrix.gz"
            
            computeMatrix reference-point \
                -S ${BW_FILES} \
                -R "${FILTERED_PEAKS}" \
                --referencePoint center \
                -a 3000 -b 3000 \
                --skipZeros \
                --missingDataAsZero \
                -o "${PEAK_MATRIX}" \
                --samplesLabel ${LABEL_ARRAY[@]} \
                -p 8
            
            if [ -f "${PEAK_MATRIX}" ]; then
                # Peak heatmap
                plotHeatmap \
                    -m "${PEAK_MATRIX}" \
                    -out "${STATE_DIR}/${condition}_${state}_peaks_heatmap.pdf" \
                    --colorMap RdBu_r \
                    --whatToShow 'plot, heatmap and colorbar' \
                    --heatmapHeight 12 \
                    --heatmapWidth 8 \
                    --plotTitle "${state} Associated Peaks - ${condition} Signal" \
                    --xAxisLabel "Distance from peak center (bp)" \
                    --yAxisLabel "Peaks (n = ${PEAK_COUNT})" \
                    --refPointLabel "Peak center" \
                    --sortRegions descend \
                    --sortUsing mean \
                    --zMin 0 --zMax auto
                
                # Peak profile plot
                plotProfile \
                    -m "${PEAK_MATRIX}" \
                    -out "${STATE_DIR}/${condition}_${state}_peaks_profile.pdf" \
                    --plotTitle "${state} Associated Peaks - ${condition} Average Signal" \
                    --yAxisLabel "Normalized read density" \
                    --xAxisLabel "Distance from peak center (bp)" \
                    --refPointLabel "Peak center" \
                    --colors "#E41A1C" "#4DAF4A" "#984EA3" \
                    --plotHeight 6 \
                    --plotWidth 8 \
                    --legendLocation upper-right \
                    --perGroup
                
                echo "      ✓ Peak heatmap and profile created"
            fi
        fi
    done
    
    # =============================
    # Create combined visualization for all conditions
    # =============================
    echo "  Creating combined visualization for ${state}..."
    
    ALL_BW=""
    ALL_LABELS=""
    
    for condition in "${CONDITIONS[@]}"; do
        for hist in "${MARKERS[@]}"; do
            bw_file="${projPath}/alignment/bigwig/human/${condition}/${hist}/normalized.bw"
            if [ -f "$bw_file" ]; then
                ALL_BW="${ALL_BW} ${bw_file}"
                ALL_LABELS="${ALL_LABELS} ${condition}_${hist}"
            fi
        done
    done
    
    if [ ! -z "$ALL_BW" ]; then
        read -ra ALL_LABEL_ARRAY <<< "${ALL_LABELS}"
        
        # Combined TSS matrix for this cell state
        computeMatrix reference-point \
            -S ${ALL_BW} \
            -R "${TSS_BED}" \
            --referencePoint TSS \
            -a 3000 -b 3000 \
            --skipZeros \
            --missingDataAsZero \
            -o ${STATE_DIR}/ALL_conditions_TSS_matrix.gz \
            --samplesLabel ${ALL_LABEL_ARRAY[@]} \
            -p 8
        
        plotHeatmap \
            -m ${STATE_DIR}/ALL_conditions_TSS_matrix.gz \
            -out ${STATE_DIR}/${state}_ALL_conditions_TSS_comparison.pdf \
            --colorMap RdBu_r \
            --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 15 \
            --heatmapWidth 14 \
            --plotTitle "${state} Cell State - All Conditions & Antibodies TSS" \
            --xAxisLabel "Distance from TSS (bp)" \
            --yAxisLabel "Genes (n = ${NUM_GENES})" \
            --refPointLabel "TSS" \
            --sortUsing mean \
            --sortRegions descend \
            --zMin 0 --zMax auto
        
        plotProfile \
            -m ${STATE_DIR}/ALL_conditions_TSS_matrix.gz \
            -out ${STATE_DIR}/${state}_ALL_conditions_TSS_profile.pdf \
            --plotTitle "${state} Cell State - TSS Signal All Conditions" \
            --yAxisLabel "Normalized read density" \
            --xAxisLabel "Distance from TSS (bp)" \
            --refPointLabel "TSS" \
            --plotHeight 8 \
            --plotWidth 12 \
            --legendLocation upper-right \
            --perGroup
        
        # Combined gene body analysis
        computeMatrix scale-regions \
            -S ${ALL_BW} \
            -R "${SUBSET_BED}" \
            --beforeRegionStartLength 3000 \
            --regionBodyLength 5000 \
            --afterRegionStartLength 3000 \
            --skipZeros \
            --missingDataAsZero \
            -o ${STATE_DIR}/ALL_conditions_gene_body_matrix.gz \
            --samplesLabel ${ALL_LABEL_ARRAY[@]} \
            -p 8
        
        plotProfile \
            -m ${STATE_DIR}/ALL_conditions_gene_body_matrix.gz \
            -out ${STATE_DIR}/${state}_ALL_conditions_gene_body_profile.pdf \
            --plotTitle "${state} Cell State - Gene Body Signal All Conditions" \
            --yAxisLabel "Normalized read density" \
            --startLabel "Gene Start" \
            --endLabel "Gene End" \
            --plotHeight 8 \
            --plotWidth 12 \
            --legendLocation upper-right \
            --perGroup
        
        echo "  ✓ Combined visualizations created"
    fi
    
    # =============================
    # Cell state-specific summary
    # =============================
    SUMMARY_FILE="${STATE_DIR}/${state}_analysis_summary.txt"
    {
        echo "Cell State Condition-Antibody Analysis Summary"
        echo "================================================"
        echo "Cell State: ${state}"
        echo "Total genes analyzed: ${NUM_GENES}"
        echo ""
        echo "Conditions analyzed: ${CONDITIONS[@]}"
        echo "Antibodies compared: ${MARKERS[@]}"
        echo ""
        if [ -s "${FILTERED_PEAKS}" ]; then
            echo "Peak Analysis:"
            echo "  Union peaks near ${state} genes: ${PEAK_COUNT}"
            echo "  Peak assignment window: ${MAX_DISTANCE} bp"
        else
            echo "Peak Analysis: No peaks found near ${state} genes"
        fi
        echo ""
        echo "Files created:"
        echo "  Per condition: TSS, gene body, and peak heatmaps/profiles"
        echo "  Combined: All conditions comparison"
    } > "${SUMMARY_FILE}"
    
    echo "  ✓ Summary report created"
done

# =============================
# Cross-cell state comparison
# =============================
echo ""
echo "================================================================"
echo "Creating cross-cell state comparisons..."
echo "================================================================"

# For each condition, create a comparison across all cell states
for condition in "${CONDITIONS[@]}"; do
    echo "  Creating comparison for ${condition}..."
    
    COMP_DIR="${OUT_DIR}/cross_state_comparison"
    mkdir -p "${COMP_DIR}"
    
    # Combine TSS beds from all cell states
    COMBINED_TSS="${COMP_DIR}/${condition}_all_states_tss.bed"
    > "${COMBINED_TSS}"
    
    for state in "${!GENE_LISTS[@]}"; do
        STATE_TSS="${OUT_DIR}/${state}/tss.bed"
        if [ -f "${STATE_TSS}" ]; then
            awk -v state="${state}" 'BEGIN{OFS="\t"} {$4=$4"_"state; print}' "${STATE_TSS}" >> "${COMBINED_TSS}"
        fi
    done
    
    if [ -s "${COMBINED_TSS}" ]; then
        sort -k1,1 -k2,2n "${COMBINED_TSS}" > "${COMBINED_TSS}.sorted"
        mv "${COMBINED_TSS}.sorted" "${COMBINED_TSS}"
        
        # Build BigWig list
        BW_FILES=""
        LABELS=""
        
        for hist in "${MARKERS[@]}"; do
            bw_file="${projPath}/alignment/bigwig/human/${condition}/${hist}/normalized.bw"
            if [ -f "$bw_file" ]; then
                BW_FILES="${BW_FILES} ${bw_file}"
                LABELS="${LABELS} ${hist}"
            fi
        done
        
        if [ ! -z "${BW_FILES}" ]; then
            read -ra LABEL_ARRAY <<< "${LABELS}"
            
            computeMatrix reference-point \
                -S ${BW_FILES} \
                -R "${COMBINED_TSS}" \
                --referencePoint TSS \
                -a 3000 -b 3000 \
                --skipZeros \
                --missingDataAsZero \
                -o ${COMP_DIR}/${condition}_all_states_TSS_matrix.gz \
                --samplesLabel ${LABEL_ARRAY[@]} \
                -p 8
            
            plotHeatmap \
                -m ${COMP_DIR}/${condition}_all_states_TSS_matrix.gz \
                -out ${COMP_DIR}/${condition}_all_states_TSS_comparison.pdf \
                --colorMap RdBu_r \
                --whatToShow 'plot, heatmap and colorbar' \
                --heatmapHeight 15 \
                --heatmapWidth 10 \
                --plotTitle "All Cell States - ${condition} TSS Signal" \
                --xAxisLabel "Distance from TSS (bp)" \
                --yAxisLabel "Genes (ADRN, MES, Int)" \
                --refPointLabel "TSS" \
                --sortUsing mean \
                --sortRegions descend \
                --zMin 0 --zMax auto
        fi
    fi
done

echo ""
echo "================================================================"
echo "ANALYSIS COMPLETE!"
echo "================================================================"
echo ""
echo "Results saved in: ${OUT_DIR}/"
echo ""
echo "For each cell state (ADRN, MES, Int):"
echo "  Per condition (Untreated, D7_Cis, D14_Cis, D14_C70):"
echo "    - {condition}_{state}_TSS_heatmap.pdf"
echo "    - {condition}_{state}_TSS_profile.pdf"
echo "    - {condition}_{state}_gene_body_heatmap.pdf"
echo "    - {condition}_{state}_gene_body_profile.pdf"
echo "    - {condition}_{state}_peaks_heatmap.pdf (if peaks found)"
echo "    - {condition}_{state}_peaks_profile.pdf (if peaks found)"
echo ""
echo "  Combined analysis:"
echo "    - {state}_ALL_conditions_TSS_comparison.pdf"
echo "    - {state}_ALL_conditions_TSS_profile.pdf"
echo "    - {state}_ALL_conditions_gene_body_profile.pdf"
echo ""
echo "Cross-state comparisons:"
echo "    - {condition}_all_states_TSS_comparison.pdf"
echo ""
echo "Key features:"
echo "  ✓ Each heatmap shows H3K4me3, KDM5A, KDM5B side-by-side"
echo "  ✓ Union peaks used for fair comparison across antibodies"
echo "  ✓ Cell state-specific regions analyzed"
echo "  ✓ TSS, gene body, and peak-centered views included"
echo "  ✓ Direct comparison of all markers at each condition"
echo "  ✓ Cross-cell state comparisons for each condition"
echo "  ✓ Color coding: H3K4me3 (red), KDM5A (green), KDM5B (purple)"