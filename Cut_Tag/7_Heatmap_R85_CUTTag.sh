#!/bin/bash
#SBATCH --job-name=TimeSeries_TSS_Peak
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
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

# Define paths
projPath="/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise"
REF_DIR="${projPath}/Reference"
OUT_DIR="${projPath}/heatmap/final_timeseries"

# Define distances to test
distanceList=(5000 7500 10000)

# Create output directory
mkdir -p ${OUT_DIR}

echo "================================================================"
echo "Creating Time Series Heatmaps: TSS and Peak-centered"
echo "2 plots per marker, fixed positions across time"
echo "================================================================"

# Step 1: Create TSS file from existing genes.bed
GENES_FILE="${projPath}/genes.bed"
TSS_FILE="${projPath}/genes_tss.bed"

if [ -f "$GENES_FILE" ]; then
    echo "Creating TSS file from genes.bed..."
    # Extract TSS: robust handling of strand and coordinates; ensure 6 fields
    awk 'BEGIN{OFS="\t"} {
        if(NF < 6) next;
        if($6=="+") {
            start=$2;
        } else if($6=="-") {
            start=$3-1;
        } else {
            next;
        }
        if(start<0) start=0;
        print $1,start,start+1,$4"_TSS",0,$6
    }' ${GENES_FILE} | sort -k1,1 -k2,2n | uniq > ${TSS_FILE}
    echo "TSS file created: $(wc -l < ${TSS_FILE}) regions"
else
    echo "ERROR: genes.bed file not found at ${GENES_FILE}"
    exit 1
fi

# Process each marker and distance
for hist in H3K4me3 KDM5A KDM5B; do
    
    echo ""
    echo "Processing ${hist}..."
    echo "========================================"
    
    # Test each distance
    for distance in ${distanceList[@]}; do
        
        echo ""
        echo "  Testing distance: ${distance} bp"
        echo "  ----------------------------------------"
        
        # Get merged peak file for this distance
        PEAK_FILE="${projPath}/Homer/mergePeaks/human/Untreated_control/${hist}/${distance}/mergedPeakFile.bed"
    
        if [ ! -f "$PEAK_FILE" ]; then
            echo "  ERROR: Missing Homer mergePeaks file for ${hist} at ${distance}bp"
            echo "         Expected: ${PEAK_FILE}"
            echo "         Please run 5_MergePeaks_R85.sh first. Skipping."
            continue
        fi
    
        PEAK_COUNT=$(wc -l < "$PEAK_FILE" 2>/dev/null || echo "0")
        echo "  Peak regions: ${PEAK_COUNT}"
    
    # ============================================
    # MERGED REPLICATE ANALYSIS (PRIMARY)
    # ============================================
    echo "  Creating MERGED replicate plots..."
    
    # Check for merged BigWig files
    BW_MERGED=""
    for condition in Untreated D7_Cis D14_Cis D14_C70; do
        bw_file="${projPath}/alignment/bigwig/human/${condition}/${hist}/normalized.bw"
        if [ -f "$bw_file" ]; then
            BW_MERGED="${BW_MERGED} ${bw_file}"
        else
            echo "    Warning: Missing merged BigWig for ${condition}"
        fi
    done
    
    if [ ! -z "$BW_MERGED" ]; then
        
        # TSS-CENTERED MERGED
        echo "    Computing TSS matrix (merged)..."
        computeMatrix reference-point \
            -S ${BW_MERGED} \
            -R ${TSS_FILE} \
            --referencePoint TSS \
            -a 3000 -b 3000 \
            --skipZeros \
            --missingDataAsZero \
            -o ${OUT_DIR}/${hist}_TSS_merged_${distance}bp_matrix.gz \
            --samplesLabel "Day 0" "Day 7" "Day 14" "Day 14+C70" \
            -p 8
        
        plotHeatmap \
            -m ${OUT_DIR}/${hist}_TSS_merged_${distance}bp_matrix.gz \
            -out ${OUT_DIR}/${hist}_TSS_merged_${distance}bp.pdf \
            --colorMap RdBu_r \
            --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 12 \
            --heatmapWidth 8 \
            --plotTitle "${hist} - TSS Signal (${distance}bp Merged)" \
            --xAxisLabel "Distance from TSS (bp)" \
            --yAxisLabel "Genes" \
            --refPointLabel "TSS" \
            --sortUsing mean \
            --sortRegions descend \
            --zMin 0 --zMax auto
        
        # PEAK-CENTERED MERGED
        echo "    Computing Peak matrix (merged)..."
        computeMatrix reference-point \
            -S ${BW_MERGED} \
            -R ${PEAK_FILE} \
            --referencePoint center \
            -a 3000 -b 3000 \
            --skipZeros \
            --missingDataAsZero \
            -o ${OUT_DIR}/${hist}_Peak_merged_${distance}bp_matrix.gz \
            --samplesLabel "Day 0" "Day 7" "Day 14" "Day 14+C70" \
            -p 8
        
        plotHeatmap \
            -m ${OUT_DIR}/${hist}_Peak_merged_${distance}bp_matrix.gz \
            -out ${OUT_DIR}/${hist}_Peak_merged_${distance}bp.pdf \
            --colorMap RdBu_r \
            --whatToShow 'plot, heatmap and colorbar' \
            --heatmapHeight 12 \
            --heatmapWidth 8 \
            --plotTitle "${hist} - Peak Signal (${distance}bp Merged)" \
            --xAxisLabel "Distance from Peak Center (bp)" \
            --yAxisLabel "Peaks (n=${PEAK_COUNT})" \
            --refPointLabel "Center" \
            --sortUsing mean \
            --sortRegions descend \
            --zMin 0 --zMax auto
    fi
    
        # ============================================
        # SEPARATE REPLICATE ANALYSIS
        # ============================================
        echo "    Creating SEPARATE replicate plots..."
        
        for rep in 1 2; do
            echo "      Replicate ${rep}..."
            
            # Check for replicate BigWig files
            BW_REP=""
            for condition in Untreated D7_Cis D14_Cis D14_C70; do
                bw_file="${projPath}/alignment/bigwig/human/${condition}/${hist}/normalized_rep${rep}.bw"
                if [ -f "$bw_file" ]; then
                    BW_REP="${BW_REP} ${bw_file}"
                fi
            done
            
            if [ ! -z "$BW_REP" ]; then
                
                # TSS-CENTERED REP
                computeMatrix reference-point \
                    -S ${BW_REP} \
                    -R ${TSS_FILE} \
                    --referencePoint TSS \
                    -a 3000 -b 3000 \
                    --skipZeros \
                    --missingDataAsZero \
                    -o ${OUT_DIR}/${hist}_TSS_rep${rep}_${distance}bp_matrix.gz \
                    --samplesLabel "Day 0" "Day 7" "Day 14" "Day 14+C70" \
                    -p 8
                
                plotHeatmap \
                    -m ${OUT_DIR}/${hist}_TSS_rep${rep}_${distance}bp_matrix.gz \
                    -out ${OUT_DIR}/${hist}_TSS_rep${rep}_${distance}bp.pdf \
                    --colorMap RdBu_r \
                    --whatToShow 'plot, heatmap and colorbar' \
                    --heatmapHeight 10 \
                    --heatmapWidth 6 \
                    --plotTitle "${hist} Rep${rep} - TSS (${distance}bp)" \
                    --xAxisLabel "Distance from TSS (bp)" \
                    --yAxisLabel "Genes" \
                    --refPointLabel "TSS" \
                    --sortUsing mean \
                    --sortRegions descend \
                    --zMin 0 --zMax auto
                
                # PEAK-CENTERED REP
                computeMatrix reference-point \
                    -S ${BW_REP} \
                    -R ${PEAK_FILE} \
                    --referencePoint center \
                    -a 3000 -b 3000 \
                    --skipZeros \
                    --missingDataAsZero \
                    -o ${OUT_DIR}/${hist}_Peak_rep${rep}_${distance}bp_matrix.gz \
                    --samplesLabel "Day 0" "Day 7" "Day 14" "Day 14+C70" \
                    -p 8
                
                plotHeatmap \
                    -m ${OUT_DIR}/${hist}_Peak_rep${rep}_${distance}bp_matrix.gz \
                    -out ${OUT_DIR}/${hist}_Peak_rep${rep}_${distance}bp.pdf \
                    --colorMap RdBu_r \
                    --whatToShow 'plot, heatmap and colorbar' \
                    --heatmapHeight 10 \
                    --heatmapWidth 6 \
                    --plotTitle "${hist} Rep${rep} - Peaks (${distance}bp)" \
                    --xAxisLabel "Distance from Peak Center (bp)" \
                    --yAxisLabel "Peaks (n=${PEAK_COUNT})" \
                    --refPointLabel "Center" \
                    --sortUsing mean \
                    --sortRegions descend \
                    --zMin 0 --zMax auto
            fi
        done
        
        echo "    ✓ Distance ${distance}bp completed for ${hist}"
    done
    
    echo "  ✓ All distances completed for ${hist}"
done

echo ""
echo "================================================================"
echo "MULTI-DISTANCE ANALYSIS COMPLETE!"
echo "================================================================"
echo ""
echo "Results in: ${OUT_DIR}/"
echo ""
echo "Files created per marker and distance:"
echo ""
for hist in H3K4me3 KDM5A KDM5B; do
    echo "  ${hist}:"
    for distance in ${distanceList[@]}; do
        echo "    Distance ${distance}bp:"
        echo "      MERGED REPLICATES:"
        echo "        - ${hist}_TSS_merged_${distance}bp.pdf    : TSS-centered, all time points"
        echo "        - ${hist}_Peak_merged_${distance}bp.pdf   : Peak-centered, all time points"
        echo "      SEPARATE REPLICATES:"
        echo "        - ${hist}_TSS_rep1_${distance}bp.pdf      : TSS-centered, replicate 1"
        echo "        - ${hist}_TSS_rep2_${distance}bp.pdf      : TSS-centered, replicate 2"
        echo "        - ${hist}_Peak_rep1_${distance}bp.pdf     : Peak-centered, replicate 1"
        echo "        - ${hist}_Peak_rep2_${distance}bp.pdf     : Peak-centered, replicate 2"
    done
done
echo ""
echo "Key features:"
echo "  ✓ Each row = same genomic position across all time points"
echo "  ✓ Columns show progression: Day 0 → Day 7 → Day 14 → Day 14+C70"
echo "  ✓ TSS plots: Signal at gene transcription start sites"
echo "  ✓ Peak plots: Signal at binding/modification sites"
echo "  ✓ Signal plot (top): Average signal across all positions"
echo "  ✓ Heatmap (bottom): Individual position signals"
echo "  ✓ Color intensity shows signal strength changes over time"
echo "  ✓ Multiple distances tested: 5000bp, 7500bp, 10000bp"
echo "  ✓ Compare distances to find optimal peak merging parameters"
echo ""
echo "Distance comparison guide:"
echo "  - 5000bp: Sharpest signals, more individual peaks"
echo "  - 7500bp: Balanced approach (recommended)"
echo "  - 10000bp: Broader regions, fewer merged peaks"
echo ""
echo "Next steps:"
echo "  1. Compare peak counts across distances"
echo "  2. Assess signal sharpness and coverage"
echo "  3. Check reproducibility across replicates"
echo "  4. Choose optimal distance for final analysis"