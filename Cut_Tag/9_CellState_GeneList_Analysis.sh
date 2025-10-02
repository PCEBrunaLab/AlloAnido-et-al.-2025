#!/bin/bash
#SBATCH --job-name=CellState_Peak_Analysis
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
OUT_DIR="${projPath}/heatmap/cellstate_genelist"
GENES_BED="${projPath}/genes.bed"

# Cell state gene files
ADRN_GENE_FILE="${projPath}/Reference/vanGroningen_2017_ADRN.csv"
MES_GENE_FILE="${projPath}/Reference/vanGroningen_2017_MES.csv"
INT_GENE_FILE="${projPath}/Reference/Int.csv"

# Distance for Homer merged peaks
distanceList=(7500)

# Peak assignment method
PEAK_ASSIGNMENT_METHOD="nearest_gene"

# Maximum distance for peak-to-gene assignment (5kb for all markers)
declare -A MAX_DISTANCE
MAX_DISTANCE["H3K4me3"]=5000
MAX_DISTANCE["KDM5A"]=5000
MAX_DISTANCE["KDM5B"]=5000
DEFAULT_MAX_DISTANCE=5000

# Markers to include
markers=(H3K4me3 KDM5A KDM5B)

# Time points
TIME_POINTS=("Untreated" "D7_Cis" "D14_Cis" "D14_C70")
TIME_LABELS=("Untreated" "Day_7" "Day_14" "Day_14+C70")

# =============================
# Setup
# =============================
mkdir -p "${OUT_DIR}"

echo "================================================================"
echo "Cell State Gene-Specific Peak and TSS Analysis"
echo "Method: ${PEAK_ASSIGNMENT_METHOD}"
echo "Peak-to-gene assignment window: 5 kb"
echo "Cell states: ADRN, MES, Int"
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
    
    echo "Debug: Creating gene subset for ${state}"
    echo "Debug: Gene list file: ${GENE_LIST} ($(wc -l < "${GENE_LIST}") genes)"
    echo "Debug: Genes annotation file: ${GENES_BED} ($(wc -l < "${GENES_BED}") total genes)"
    
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]; next} ($4 in a)' \
        "${GENE_LIST}" "${GENES_BED}" | \
        sort -k1,1 -k2,2n | uniq > "${SUBSET_BED}"
    
    NUM_GENES=$(wc -l < "${SUBSET_BED}" 2>/dev/null || echo "0")
    if [ "${NUM_GENES}" -eq 0 ]; then
        echo "WARNING: No matching genes found for ${state}. Skipping."
        echo "Debug: Check if gene names in ${GENE_LIST} match gene names in ${GENES_BED}"
        echo "Debug: Sample genes from list: $(head -3 "${GENE_LIST}" | tr '\n' ' ')"
        echo "Debug: Sample genes from annotation: $(head -3 "${GENES_BED}" | cut -f4 | tr '\n' ' ')"
        continue
    fi
    
    echo "Matched genes: ${NUM_GENES}"
    
    # Create TSS BED file
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
    
    # Process each marker
    for hist in "${markers[@]}"; do
        echo ""
        echo "Processing ${hist} for ${state}"
        echo "========================================"
        
        # Build BigWig file lists
        BW_FILES=""
        ACTUAL_LABELS=""
        
        for i in "${!TIME_POINTS[@]}"; do
            tp="${TIME_POINTS[$i]}"
            label="${TIME_LABELS[$i]}"
            bw_file="${projPath}/alignment/bigwig/human/${tp}/${hist}/normalized.bw"
            
            if [ -f "$bw_file" ]; then
                BW_FILES="${BW_FILES} ${bw_file}"
                ACTUAL_LABELS="${ACTUAL_LABELS} ${label}"
            else
                echo "    Warning: Missing BigWig for ${tp}"
            fi
        done
        
        if [ -z "${BW_FILES}" ]; then
            echo "    ERROR: No BigWig files found for ${hist}"
            continue
        fi
        
        read -ra LABEL_ARRAY <<< "${ACTUAL_LABELS}"
        
        # TSS-centered analysis
        echo "  Creating TSS-centered heatmap..."
        TSS_MATRIX="${STATE_DIR}/${hist}_TSS_matrix.gz"
        
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
            plotHeatmap \
                -m "${TSS_MATRIX}" \
                -out "${STATE_DIR}/${hist}_${state}_TSS_heatmap.pdf" \
                --colorMap RdBu_r \
                --whatToShow 'plot, heatmap and colorbar' \
                --heatmapHeight 15 \
                --heatmapWidth 10 \
                --plotTitle "${hist} Signal at ${state} Cell State Gene TSS" \
                --xAxisLabel "Distance from TSS (bp)" \
                --yAxisLabel "Genes (n = ${NUM_GENES})" \
                --refPointLabel "TSS" \
                --sortRegions descend \
                --sortUsing mean \
                --zMin 0 --zMax auto
            
            plotProfile \
                -m "${TSS_MATRIX}" \
                -out "${STATE_DIR}/${hist}_${state}_TSS_profile.pdf" \
                --plotTitle "${hist} Average Signal at ${state} Cell State Gene TSS" \
                --yAxisLabel "Normalized read density" \
                --refPointLabel "TSS" \
                --colors "#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" \
                --plotHeight 8 \
                --plotWidth 10 \
                --legendLocation upper-right \
                --perGroup
            
            echo "    ✓ TSS heatmap and profile created"
        fi
        
        # =============================
        # Gene Body analysis (NEW - captures entire gene body)
        # =============================
        echo "  Creating GENE BODY-scaled heatmaps..."
        
        # Use scale-regions instead of reference-point
        # This captures the ENTIRE GENE from start to end
        GENE_BODY_MATRIX="${STATE_DIR}/${hist}_gene_body_matrix.gz"
        
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
                -out "${STATE_DIR}/${hist}_${state}_gene_body_heatmap.pdf" \
                --colorMap RdBu_r \
                --whatToShow 'plot, heatmap and colorbar' \
                --heatmapHeight 15 \
                --heatmapWidth 10 \
                --plotTitle "${hist} Signal Across ${state} Gene Bodies" \
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
                -out "${STATE_DIR}/${hist}_${state}_gene_body_profile.pdf" \
                --plotTitle "${hist} Average Signal Across ${state} Gene Bodies" \
                --yAxisLabel "Normalized read density" \
                --startLabel "Gene Start" \
                --endLabel "Gene End" \
                --colors "#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" \
                --plotHeight 8 \
                --plotWidth 10 \
                --legendLocation upper-right \
                --perGroup
            
            echo "    ✓ Gene body heatmap and profile created"
        fi
        
        # Peak-centered analysis
        for distance in "${distanceList[@]}"; do
            echo "  Homer merge distance: ${distance} bp"
            
            PEAK_FILE="${projPath}/Homer/mergePeaks/human/Untreated_control/${hist}/${distance}/mergedPeakFile.bed"
            
            if [ ! -f "${PEAK_FILE}" ]; then
                echo "    Peak file not found: ${PEAK_FILE}"
                continue
            fi
            
            # Check chromosome format consistency
            GENE_CHR=$(head -1 "${SUBSET_BED}" 2>/dev/null | cut -f1 | grep -c "^chr" || echo "0")
            PEAK_CHR=$(head -1 "${PEAK_FILE}" 2>/dev/null | cut -f1 | grep -c "^chr" || echo "0")

            PEAK_CLEAN="${STATE_DIR}/${hist}_peaks_clean_${distance}bp.bed"

            # Debug: Check formats
            echo "    Debug: Gene chr format (1=has chr): ${GENE_CHR}"
            echo "    Debug: Peak chr format (1=has chr): ${PEAK_CHR}"

            if [ "${GENE_CHR}" != "${PEAK_CHR}" ]; then
                echo "    Fixing chromosome naming mismatch..."
                if [ "${GENE_CHR}" -eq 1 ]; then
                    # Genes have chr, add to peaks if missing
                    awk 'BEGIN{OFS="\t"} {
                        if(NF<3) next;
                        if($2<0) $2=0;
                        if($3<=$2) $3=$2+1;
                        if($1 !~ /^chr/) $1="chr"$1;
                        print $1,$2,$3,(NF>=4 ? $4 : "peak_"NR)
                    }' "${PEAK_FILE}" | sort -k1,1 -k2,2n > "${PEAK_CLEAN}"
                else
                    # Genes don't have chr, remove from peaks if present
                    awk 'BEGIN{OFS="\t"} {
                        if(NF<3) next;
                        if($2<0) $2=0;
                        if($3<=$2) $3=$2+1;
                        gsub(/^chr/, "", $1);
                        print $1,$2,$3,(NF>=4 ? $4 : "peak_"NR)
                    }' "${PEAK_FILE}" | sort -k1,1 -k2,2n > "${PEAK_CLEAN}"
                fi
            else
                # No mismatch, just clean the file
                awk 'BEGIN{OFS="\t"} {
                    if(NF<3) next;
                    if($2<0) $2=0;
                    if($3<=$2) $3=$2+1;
                    print $1,$2,$3,(NF>=4 ? $4 : "peak_"NR)
                }' "${PEAK_FILE}" | sort -k1,1 -k2,2n > "${PEAK_CLEAN}"
            fi

            # Debug: Check if files exist and have content
            echo "    Debug: Peak file lines: $(wc -l < "${PEAK_FILE}" 2>/dev/null || echo "0")"
            echo "    Debug: Clean peak lines: $(wc -l < "${PEAK_CLEAN}" 2>/dev/null || echo "0")"
            echo "    Debug: Subset BED lines: $(wc -l < "${SUBSET_BED}" 2>/dev/null || echo "0")"

            # Debug: Show sample of each file
            echo "    Debug: First peak:"
            head -1 "${PEAK_CLEAN}" 2>/dev/null || echo "    No peaks"
            echo "    Debug: First gene:"
            head -1 "${SUBSET_BED}" 2>/dev/null || echo "    No genes"
            
            # Gene-Specific Peak Assignment
            echo "    Assigning peaks to nearest genes within 5 kb"
            
            MAX_DIST="${MAX_DISTANCE[$hist]:-$DEFAULT_MAX_DISTANCE}"
            
            FILTERED_PEAKS="${STATE_DIR}/${hist}_gene_associated_peaks_${distance}.bed"
            
            bedtools window -w ${MAX_DIST} -a "${PEAK_CLEAN}" -b "${SUBSET_BED}" | \
                awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_near_"$8}' | \
                sort -u > "${FILTERED_PEAKS}"
            
            # Debug: Check if peaks were found
            echo "    Debug: Filtered peaks found: $(wc -l < "${FILTERED_PEAKS}" 2>/dev/null || echo "0")"
            
            # If no peaks found, do a test with larger window
            if [ ! -s "${FILTERED_PEAKS}" ]; then
                echo "    Testing with 50kb window..."
                TEST_COUNT=$(bedtools window -w 50000 -a "${PEAK_CLEAN}" -b "${SUBSET_BED}" | wc -l)
                echo "    Peaks within 50kb: ${TEST_COUNT}"
                
                # Show sample of what bedtools window returns
                echo "    Debug: Sample bedtools window output (50kb):"
                bedtools window -w 50000 -a "${PEAK_CLEAN}" -b "${SUBSET_BED}" | head -3 2>/dev/null || echo "    No overlaps found"
                
                echo "    No peaks found within ${MAX_DIST} bp of ${state} genes"
                continue
            fi
            
            PEAK_COUNT=$(wc -l < "${FILTERED_PEAKS}")
            echo "    Found ${PEAK_COUNT} peaks associated with ${state} genes"
            
            # Compute matrix for peaks
            PEAK_MATRIX="${STATE_DIR}/${hist}_peaks_matrix.gz"
            
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
                # Plot heatmap with sorting
                plotHeatmap \
                    -m "${PEAK_MATRIX}" \
                    -out "${STATE_DIR}/${hist}_${state}_peaks_heatmap.pdf" \
                    --colorMap RdBu_r \
                    --whatToShow 'plot, heatmap and colorbar' \
                    --heatmapHeight 15 \
                    --heatmapWidth 10 \
                    --plotTitle "${hist} Signal at ${state}-Associated Peaks" \
                    --xAxisLabel "Distance from peak center (bp)" \
                    --yAxisLabel "Genes (n = ${PEAK_COUNT})" \
                    --refPointLabel "Peak center" \
                    --sortRegions descend \
                    --sortUsing mean \
                    --zMin 0 --zMax auto
                
                # Profile plot
                plotProfile \
                    -m "${PEAK_MATRIX}" \
                    -out "${STATE_DIR}/${hist}_${state}_peaks_profile.pdf" \
                    --plotTitle "${hist} Average Signal at ${state}-Associated Peaks" \
                    --yAxisLabel "Normalized read density" \
                    --refPointLabel "Peak center" \
                    --colors "#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" \
                    --plotHeight 8 \
                    --plotWidth 10 \
                    --legendLocation upper-right \
                    --perGroup
                
                echo "    ✓ Peak heatmap and profile created"
            fi
            
            # Create assignment summary
            SUMMARY_FILE="${STATE_DIR}/${hist}_peak_assignment_summary.txt"
            echo "Peak-Gene Assignment Summary" > "${SUMMARY_FILE}"
            echo "=============================" >> "${SUMMARY_FILE}"
            echo "Cell State: ${state}" >> "${SUMMARY_FILE}"
            echo "Chromatin mark: ${hist}" >> "${SUMMARY_FILE}"
            echo "Assignment window: ${MAX_DIST} bp" >> "${SUMMARY_FILE}"
            echo "Total peaks assigned: ${PEAK_COUNT}" >> "${SUMMARY_FILE}"
            echo "" >> "${SUMMARY_FILE}"
            
            # Distance distribution
            echo "Distance distribution:" >> "${SUMMARY_FILE}"
            awk -F'_near_' '{print $2}' "${FILTERED_PEAKS}" | \
                awk -F'_dist' '{print $2}' | \
                awk '{if($1<=2000) p++; else if($1<=5000) d++} 
                     END {print "  Promoter (≤2kb): " p "\n  Distal (2-5kb): " d}' >> "${SUMMARY_FILE}"
            
            # Top genes with most peaks
            echo "" >> "${SUMMARY_FILE}"
            echo "Genes with multiple peaks:" >> "${SUMMARY_FILE}"
            awk -F'_near_' '{print $2}' "${FILTERED_PEAKS}" | \
                awk -F'_dist' '{print $1}' | \
                sort | uniq -c | sort -rn | head -10 | \
                awk '{print "  " $2 ": " $1 " peaks"}' >> "${SUMMARY_FILE}"
            
            echo "    ✓ Summary created"
        done
    done
done

# =============================
# ACTIVE GENE ANALYSIS SECTION
# =============================
echo ""
echo "================================================================"
echo "EXTRACTING ACTIVE GENE SETS BASED ON H3K4ME3"
echo "================================================================"

for state in "${!GENE_LISTS[@]}"; do
    STATE_DIR="${OUT_DIR}/${state}"
    
    if [ ! -f "${STATE_DIR}/H3K4me3_gene_associated_peaks_7500.bed" ]; then
        echo "Skipping ${state} - no H3K4me3 peak data"
        continue
    fi
    
    echo ""
    echo "Processing active genes for ${state} cell state..."
    
    # Extract genes with H3K4me3 peaks
    ACTIVE_GENES="${STATE_DIR}/H3K4me3_active_genes.txt"
    awk -F'_near_' '{print $2}' "${STATE_DIR}/H3K4me3_gene_associated_peaks_7500.bed" | \
        awk -F'_dist' '{print $1}' | sort -u > "${ACTIVE_GENES}"
    
    ACTIVE_COUNT=$(wc -l < "${ACTIVE_GENES}")
    ORIGINAL_COUNT=$(wc -l < "${GENE_LISTS[${state}]}")
    
    if [ "${ORIGINAL_COUNT}" -gt 0 ]; then
        PERCENT_ACTIVE=$(echo "scale=1; ${ACTIVE_COUNT}*100/${ORIGINAL_COUNT}" | bc)
    else
        PERCENT_ACTIVE="0"
    fi
    
    echo "  Original genes: ${ORIGINAL_COUNT}"
    echo "  Active genes (with H3K4me3): ${ACTIVE_COUNT} (${PERCENT_ACTIVE}%)"
    
    # Create active gene subset BED
    SUBSET_BED="${STATE_DIR}/genes_subset.bed"
    ACTIVE_GENES_BED="${STATE_DIR}/active_genes_only.bed"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]; next} ($4 in a)' \
        "${ACTIVE_GENES}" "${SUBSET_BED}" > "${ACTIVE_GENES_BED}"
    
    # Check KDM5A/B overlap with active genes
    echo "  Checking KDM5 presence at active genes:"
    for hist in KDM5A KDM5B; do
        if [ -f "${STATE_DIR}/${hist}_gene_associated_peaks_7500.bed" ]; then
            # Get KDM5 peaks at active genes
            awk -F'_near_' 'NR==FNR{a[$1]; next} ($2 in a){print}' \
                "${ACTIVE_GENES}" \
                "${STATE_DIR}/${hist}_gene_associated_peaks_7500.bed" > \
                "${STATE_DIR}/${hist}_at_active_genes.bed"
            
            KDM5_AT_ACTIVE=$(awk -F'_near_' '{print $2}' "${STATE_DIR}/${hist}_at_active_genes.bed" | \
                            awk -F'_dist' '{print $1}' | sort -u | wc -l)
            
            if [ "${ACTIVE_COUNT}" -gt 0 ]; then
                PERCENT_WITH_KDM5=$(echo "scale=1; ${KDM5_AT_ACTIVE}*100/${ACTIVE_COUNT}" | bc)
            else
                PERCENT_WITH_KDM5="0"
            fi
            echo "    ${hist}: ${KDM5_AT_ACTIVE}/${ACTIVE_COUNT} active genes (${PERCENT_WITH_KDM5}%)"
        fi
    done
    
    # Save summary
    ACTIVE_SUMMARY="${STATE_DIR}/active_gene_summary.txt"
    {
        echo "Active Gene Analysis Summary"
        echo "============================"
        echo "Cell State: ${state}"
        echo "Total genes in list: ${ORIGINAL_COUNT}"
        echo "Active genes (H3K4me3+): ${ACTIVE_COUNT} (${PERCENT_ACTIVE}%)"
        echo ""
        echo "Top 10 active genes:"
        head -10 "${ACTIVE_GENES}"
    } > "${ACTIVE_SUMMARY}"
done

echo ""
echo "Active gene analysis complete!"
echo "Check active_gene_summary.txt in each condition folder"

echo ""
echo "================================================================"
echo "ANALYSIS COMPLETE"
echo "================================================================"
echo "Cell state analysis includes TSS, gene body, and peak-centered views"
echo "Peak assignment: Nearest gene within 5 kb"
echo "Results saved in: ${OUT_DIR}"
echo ""
echo "Files created per cell state/marker:"
echo "  - TSS heatmap and profile (±3kb around TSS)"
echo "  - Gene body heatmap and profile (full gene + flanking)"
echo "  - Peak heatmap and profile (±3kb around peaks)"
echo "  - Analysis summary"
echo "  - Active gene analysis (H3K4me3-based)"
