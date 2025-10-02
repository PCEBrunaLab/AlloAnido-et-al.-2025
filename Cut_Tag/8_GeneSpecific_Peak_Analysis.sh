#!/bin/bash
#SBATCH --job-name=GeneSpecific_Peak_Analysis
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
OUT_DIR="${projPath}/heatmap/gene_specific_peaks"
CSV_FILE="${projPath}/Reference/gene_list.csv"
GENES_BED="${projPath}/genes.bed"

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
echo "Gene-Specific Peak and TSS Analysis"
echo "Method: ${PEAK_ASSIGNMENT_METHOD}"
echo "Peak-to-gene assignment window: 5 kb"
echo "================================================================"

# Validate inputs
if [ ! -f "${GENES_BED}" ]; then
    echo "ERROR: genes.bed file not found at ${GENES_BED}"
    exit 1
fi

if [ ! -f "${CSV_FILE}" ]; then
    echo "ERROR: CSV file not found at ${CSV_FILE}"
    exit 1
fi

# =============================
# Parse CSV and create gene lists
# =============================
echo "Parsing CSV file: ${CSV_FILE}"

# Remove BOM and extract header
HEADER=$(head -n1 "${CSV_FILE}" | sed 's/^\xEF\xBB\xBF//' | tr -d '\r')
IFS=',' read -ra CONDITIONS <<< "${HEADER}"

# Clean condition names
for i in "${!CONDITIONS[@]}"; do
    CONDITIONS[i]=$(echo "${CONDITIONS[i]}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/^"//;s/"$//')
done

echo "Detected conditions: ${CONDITIONS[*]}"

# Create gene lists for each condition
declare -A GENE_LISTS
for condition in "${CONDITIONS[@]}"; do
    if [ -z "${condition}" ]; then continue; fi
    
    echo "Processing condition: ${condition}"
    
    # Find column index
    col_index=1
    for i in "${!CONDITIONS[@]}"; do
        if [ "${CONDITIONS[i]}" = "${condition}" ]; then
            col_index=$((i+1))
            break
        fi
    done
    
    # Extract genes for this condition
    GENE_LIST_FILE="${OUT_DIR}/${condition}_genes.txt"
    tail -n +2 "${CSV_FILE}" | cut -d',' -f"${col_index}" | \
        sed 's/^[[:space:]]*//;s/[[:space:]]*$//;s/^"//;s/"$//' | \
        grep -v -E '^\s*$' | sort -u > "${GENE_LIST_FILE}"
    
    NUM_GENES=$(wc -l < "${GENE_LIST_FILE}" 2>/dev/null || echo "0")
    echo "  ${condition}: ${NUM_GENES} genes"
    
    if [ "${NUM_GENES}" -gt 0 ]; then
        GENE_LISTS[${condition}]="${GENE_LIST_FILE}"
    fi
done

# =============================
# Process each condition
# =============================
for condition in "${!GENE_LISTS[@]}"; do
    echo ""
    echo "================================================================"
    echo "Processing condition: ${condition}"
    echo "================================================================"
    
    GENE_LIST="${GENE_LISTS[${condition}]}"
    CONDITION_DIR="${OUT_DIR}/${condition}"
    mkdir -p "${CONDITION_DIR}"
    
    # Create gene subset BED
    SUBSET_BED="${CONDITION_DIR}/genes_subset.bed"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]; next} ($4 in a)' \
        "${GENE_LIST}" "${GENES_BED}" | \
        sort -k1,1 -k2,2n | uniq > "${SUBSET_BED}"
    
    NUM_GENES=$(wc -l < "${SUBSET_BED}" 2>/dev/null || echo "0")
    if [ "${NUM_GENES}" -eq 0 ]; then
        echo "WARNING: No matching genes found for ${condition}. Skipping."
        continue
    fi
    
    echo "Matched genes: ${NUM_GENES}"
    
    # =============================
    # Create TSS BED file for this condition
    # =============================
    TSS_BED="${CONDITION_DIR}/tss.bed"
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
    # Process each marker
    # =============================
    for hist in "${markers[@]}"; do
        echo ""
        echo "Processing ${hist} for ${condition}"
        echo "========================================"
        
        # Build BigWig file lists (used for both TSS and peaks)
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
        
        # Convert space-separated labels to array
        read -ra LABEL_ARRAY <<< "${ACTUAL_LABELS}"
        
        # =============================
        # TSS-centered analysis
        # =============================
        echo "  Creating TSS-centered heatmap..."
        TSS_MATRIX="${CONDITION_DIR}/${hist}_TSS_matrix.gz"
        
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
                -out "${CONDITION_DIR}/${hist}_${condition}_TSS_heatmap.pdf" \
                --colorMap RdBu_r \
                --whatToShow 'plot, heatmap and colorbar' \
                --heatmapHeight 15 \
                --heatmapWidth 10 \
                --plotTitle "${hist} Signal at ${condition} Gene TSS" \
                --xAxisLabel "Distance from TSS (bp)" \
                --yAxisLabel "Genes (n = ${NUM_GENES})" \
                --refPointLabel "TSS" \
                --sortRegions descend \
                --sortUsing mean \
                --zMin 0 --zMax auto
            
            # TSS profile plot
            plotProfile \
                -m "${TSS_MATRIX}" \
                -out "${CONDITION_DIR}/${hist}_${condition}_TSS_profile.pdf" \
                --plotTitle "${hist} Average Signal at ${condition} Gene TSS" \
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
        # Peak-centered analysis
        # =============================
        for distance in "${distanceList[@]}"; do
            echo "  Homer merge distance: ${distance} bp"
            
            # Get peak file
            PEAK_FILE="${projPath}/Homer/mergePeaks/human/Untreated_control/${hist}/${distance}/mergedPeakFile.bed"
            
            if [ ! -f "${PEAK_FILE}" ]; then
                echo "    Peak file not found: ${PEAK_FILE}"
                continue
            fi
            
            # Check chromosome format consistency
            GENE_CHR=$(head -1 "${SUBSET_BED}" 2>/dev/null | cut -f1 | grep -c "^chr" || echo "0")
            PEAK_CHR=$(head -1 "${PEAK_FILE}" 2>/dev/null | cut -f1 | grep -c "^chr" || echo "0")
            
            PEAK_CLEAN="${CONDITION_DIR}/${hist}_peaks_clean_${distance}bp.bed"
            
            if [ "${GENE_CHR}" != "${PEAK_CHR}" ]; then
                echo "    Fixing chromosome naming mismatch..."
                if [ "${GENE_CHR}" -eq 1 ]; then
                    # Genes have chr, add to peaks if missing
                    awk 'BEGIN{OFS="\t"} {
                        if(NF<3) next;
                        if($2<0) $2=0;
                        if($3<=$2) $3=$2+1;
                        if($1 !~ /^chr/) $1 = "chr" $1;
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
            
            # Gene-Specific Peak Assignment using bedtools window
            echo "    Assigning peaks to nearest genes within 5 kb"
            
            # Get the 5kb distance for this marker
            MAX_DIST="${MAX_DISTANCE[$hist]:-$DEFAULT_MAX_DISTANCE}"
            
            FILTERED_PEAKS="${CONDITION_DIR}/${hist}_gene_associated_peaks_${distance}.bed"
            
            # Use bedtools window to find peaks within 5kb of our specific genes
            bedtools window -w ${MAX_DIST} -a "${PEAK_CLEAN}" -b "${SUBSET_BED}" | \
                awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4"_near_"$8}' | \
                sort -u > "${FILTERED_PEAKS}"
            
            # Check if we found any peaks
            if [ ! -s "${FILTERED_PEAKS}" ]; then
                echo "    No peaks found within ${MAX_DIST} bp of ${condition} genes"
                continue
            fi
            
            PEAK_COUNT=$(wc -l < "${FILTERED_PEAKS}")
            echo "    Found ${PEAK_COUNT} peaks associated with ${condition} genes"
            
            # Create peak-centered heatmaps
            echo "    Creating peak-centered heatmaps and profiles..."
            
            # Compute matrix
            PEAK_MATRIX="${CONDITION_DIR}/${hist}_peaks_matrix.gz"
            
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
            
            if [ ! -f "${PEAK_MATRIX}" ]; then
                echo "    ERROR: Matrix computation failed"
                continue
            fi
            
            # Plot heatmap with sorting
            plotHeatmap \
                -m "${PEAK_MATRIX}" \
                -out "${CONDITION_DIR}/${hist}_${condition}_peaks_heatmap.pdf" \
                --colorMap RdBu_r \
                --whatToShow 'plot, heatmap and colorbar' \
                --heatmapHeight 15 \
                --heatmapWidth 10 \
                --plotTitle "${hist} Signal at ${condition}-Associated Peaks" \
                --xAxisLabel "Distance from peak center (bp)" \
                --yAxisLabel "Peaks (n = ${PEAK_COUNT})" \
                --refPointLabel "Peak center" \
                --sortRegions descend \
                --sortUsing mean \
                --zMin 0 --zMax auto
            
            # Profile plot
            plotProfile \
                -m "${PEAK_MATRIX}" \
                -out "${CONDITION_DIR}/${hist}_${condition}_peaks_profile.pdf" \
                --plotTitle "${hist} Average Signal at ${condition}-Associated Peaks" \
                --yAxisLabel "Normalized read density" \
                --refPointLabel "Peak center" \
                --colors "#1f77b4" "#ff7f0e" "#2ca02c" "#d62728" \
                --plotHeight 8 \
                --plotWidth 10 \
                --legendLocation upper-right \
                --perGroup
            
            echo "    ✓ Peak heatmap and profile created"
            
            # Create assignment summary
            SUMMARY_FILE="${CONDITION_DIR}/${hist}_analysis_summary.txt"
            echo "Analysis Summary" > "${SUMMARY_FILE}"
            echo "================" >> "${SUMMARY_FILE}"
            echo "Condition: ${condition}" >> "${SUMMARY_FILE}"
            echo "Chromatin mark: ${hist}" >> "${SUMMARY_FILE}"
            echo "" >> "${SUMMARY_FILE}"
            echo "TSS Analysis:" >> "${SUMMARY_FILE}"
            echo "  Genes analyzed: ${NUM_GENES}" >> "${SUMMARY_FILE}"
            echo "" >> "${SUMMARY_FILE}"
            echo "Peak Analysis:" >> "${SUMMARY_FILE}"
            echo "  Assignment window: ${MAX_DIST} bp" >> "${SUMMARY_FILE}"
            echo "  Total peaks assigned: ${PEAK_COUNT}" >> "${SUMMARY_FILE}"
            echo "" >> "${SUMMARY_FILE}"
            
            # Top genes with most peaks
            echo "Genes with multiple peaks:" >> "${SUMMARY_FILE}"
            awk -F'_near_' '{print $2}' "${FILTERED_PEAKS}" | \
                sort | uniq -c | sort -rn | head -10 | \
                awk '{print "  " $2 ": " $1 " peaks"}' >> "${SUMMARY_FILE}"
            
            echo "    ✓ Summary created"
        done
    done
done

# =============================
# ACTIVE GENE ANALYSIS SECTION (NEW)
# =============================
echo ""
echo "================================================================"
echo "EXTRACTING ACTIVE GENE SETS BASED ON H3K4ME3"
echo "================================================================"

for condition in "${!GENE_LISTS[@]}"; do
    CONDITION_DIR="${OUT_DIR}/${condition}"
    
    # Check if H3K4me3 peak analysis was done
    if [ ! -f "${CONDITION_DIR}/H3K4me3_gene_associated_peaks_7500.bed" ]; then
        echo "Skipping ${condition} - no H3K4me3 peak data"
        continue
    fi
    
    echo ""
    echo "Processing active genes for ${condition}..."
    
    # Extract genes with H3K4me3 peaks
    ACTIVE_GENES="${CONDITION_DIR}/H3K4me3_active_genes.txt"
    awk -F'_near_' '{print $2}' "${CONDITION_DIR}/H3K4me3_gene_associated_peaks_7500.bed" | \
        sort -u > "${ACTIVE_GENES}"
    
    ACTIVE_COUNT=$(wc -l < "${ACTIVE_GENES}")
    ORIGINAL_COUNT=$(wc -l < "${GENE_LISTS[${condition}]}")
    PERCENT_ACTIVE=$(echo "scale=1; ${ACTIVE_COUNT}*100/${ORIGINAL_COUNT}" | bc)
    
    echo "  Original genes: ${ORIGINAL_COUNT}"
    echo "  Active genes (with H3K4me3): ${ACTIVE_COUNT} (${PERCENT_ACTIVE}%)"
    
    # Create active gene subset BED
    SUBSET_BED="${CONDITION_DIR}/genes_subset.bed"
    ACTIVE_GENES_BED="${CONDITION_DIR}/active_genes_only.bed"
    awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]; next} ($4 in a)' \
        "${ACTIVE_GENES}" "${SUBSET_BED}" > "${ACTIVE_GENES_BED}"
    
    # Check KDM5A/B overlap with active genes
    echo "  Checking KDM5 presence at active genes:"
    for hist in KDM5A KDM5B; do
        if [ -f "${CONDITION_DIR}/${hist}_gene_associated_peaks_7500.bed" ]; then
            # Get KDM5 peaks at active genes
            awk -F'_near_' 'NR==FNR{a[$1]; next} ($2 in a){print}' \
                "${ACTIVE_GENES}" \
                "${CONDITION_DIR}/${hist}_gene_associated_peaks_7500.bed" > \
                "${CONDITION_DIR}/${hist}_at_active_genes.bed"
            
            KDM5_AT_ACTIVE=$(cut -d'_' -f3 "${CONDITION_DIR}/${hist}_at_active_genes.bed" | \
                            sed 's/near_//' | sort -u | wc -l)
            
            PERCENT_WITH_KDM5=$(echo "scale=1; ${KDM5_AT_ACTIVE}*100/${ACTIVE_COUNT}" | bc)
            echo "    ${hist}: ${KDM5_AT_ACTIVE}/${ACTIVE_COUNT} active genes (${PERCENT_WITH_KDM5}%)"
        fi
    done
    
    # Save summary
    ACTIVE_SUMMARY="${CONDITION_DIR}/active_gene_summary.txt"
    {
        echo "Active Gene Analysis Summary"
        echo "============================"
        echo "Condition: ${condition}"
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
echo "Analysis includes both TSS and peak-centered views"
echo "Peak assignment: Nearest gene within 5 kb"
echo "Results saved in: ${OUT_DIR}"
echo ""
echo "Files created per condition/marker:"
echo "  - TSS heatmap and profile"
echo "  - Peak heatmap and profile"
echo "  - Analysis summary"
echo "  - Active gene analysis (H3K4me3-based)"