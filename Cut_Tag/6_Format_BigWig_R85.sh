#!/bin/bash
#SBATCH --job-name=BigWig_Normalize_R85
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

# Activate UCSC tools environment (has bamCoverage 3.5.6)
conda activate ucsc_tools

# Load SAMtools separately
module load SAMtools/1.11

# Verify bamCoverage is available
echo "Using bamCoverage version: $(bamCoverage --version)"

# Define paths
projPath="/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise"
bamPath="${projPath}/bam_files_dedup"
bigwigPath="${projPath}/alignment/bigwig"

# Following tutorial exactly
org="human"
gsize="2913022398"  # Human genome size for RPGC normalization
Time=("Untreated" "D7_Cis" "D14_Cis" "D14_C70")
histList=("H3K4me3" "KDM5A" "KDM5B")  # Excluding IgG as you mentioned

# Sample mapping
declare -A sample_map
# Untreated samples
sample_map["Untreated_H3K4me3_1"]="R85_0003"
sample_map["Untreated_H3K4me3_2"]="R85_0009"
sample_map["Untreated_KDM5A_1"]="R85_0015"
sample_map["Untreated_KDM5A_2"]="R85_0021"
sample_map["Untreated_KDM5B_1"]="R85_0027"
sample_map["Untreated_KDM5B_2"]="R85_0033"

# D7_Cis samples
sample_map["D7_Cis_H3K4me3_1"]="R85_0004"
sample_map["D7_Cis_H3K4me3_2"]="R85_0010"
sample_map["D7_Cis_KDM5A_1"]="R85_0016"
sample_map["D7_Cis_KDM5A_2"]="R85_0022"
sample_map["D7_Cis_KDM5B_1"]="R85_0028"
sample_map["D7_Cis_KDM5B_2"]="R85_0034"

# D14_Cis samples
sample_map["D14_Cis_H3K4me3_1"]="R85_0005"
sample_map["D14_Cis_H3K4me3_2"]="R85_0011"
sample_map["D14_Cis_KDM5A_1"]="R85_0017"
sample_map["D14_Cis_KDM5A_2"]="R85_0023"
sample_map["D14_Cis_KDM5B_1"]="R85_0029"
sample_map["D14_Cis_KDM5B_2"]="R85_0035"

# D14_C70 samples
sample_map["D14_C70_H3K4me3_1"]="R85_0006"
sample_map["D14_C70_H3K4me3_2"]="R85_0012"
sample_map["D14_C70_KDM5A_1"]="R85_0018"
sample_map["D14_C70_KDM5A_2"]="R85_0024"
sample_map["D14_C70_KDM5B_1"]="R85_0030"
sample_map["D14_C70_KDM5B_2"]="R85_0036"

echo "================================================================"
echo "Creating RPGC Normalized BigWig Files (Following Tutorial)"
echo "================================================================"

# Following tutorial loop structure exactly
for time in ${Time[@]}; do
    for hist in ${histList[@]}; do
        # Process each replicate
        for rep in 1 2; do
            sample_key="${time}_${hist}_${rep}"
            sample_id="${sample_map[$sample_key]}"
            
            if [ -z "$sample_id" ]; then
                echo "Warning: No sample mapping for ${sample_key}"
                continue
            fi
            
            echo "Processing: ${time} ${hist} rep${rep} (${sample_id})"
            
            # Create directory following tutorial structure
            mkdir -p ${bigwigPath}/${org}/${time}/${hist}
            
            # Input BAM and output BigWig
            input_bam="${bamPath}/${sample_id}_dedup.bam"
            output_bw="${bigwigPath}/${org}/${time}/${hist}/normalized_rep${rep}.bw"
            
            if [ ! -f "$input_bam" ]; then
                echo "  Warning: BAM file not found: $input_bam"
                continue
            fi
            
            # First sort and index if not already done
            sorted_bam="${bamPath}/${sample_id}_dedup.sorted.bam"
            if [ ! -f "${sorted_bam}.bai" ]; then
                echo "  Sorting and indexing BAM..."
                samtools sort -@ 8 -o ${sorted_bam} ${input_bam}
                samtools index ${sorted_bam}
            else
                sorted_bam=${input_bam}
            fi
            
            # Create normalized BigWig following tutorial exactly
            echo "  Creating RPGC normalized BigWig..."
            bamCoverage --normalizeUsing RPGC \
                --effectiveGenomeSize ${gsize} \
                -b ${sorted_bam} \
                -o ${output_bw} \
                --binSize 10 \
                --numberOfProcessors 8
            
            echo "  Created: ${output_bw}"
        done
        
        # Create merged BigWig for both replicates
        echo "  Creating merged normalized BigWig for ${time} ${hist}..."
        
        bam1="${bamPath}/${sample_map[${time}_${hist}_1]}_dedup.bam"
        bam2="${bamPath}/${sample_map[${time}_${hist}_2]}_dedup.bam"
        
        if [ -f "$bam1" ] && [ -f "$bam2" ]; then
            # Merge BAMs first
            merged_bam="${bigwigPath}/${org}/${time}/${hist}/merged.bam"
            samtools merge -@ 8 ${merged_bam} ${bam1} ${bam2}
            samtools index ${merged_bam}
            
            # Create merged normalized BigWig
            bamCoverage --normalizeUsing RPGC \
                --effectiveGenomeSize ${gsize} \
                -b ${merged_bam} \
                -o ${bigwigPath}/${org}/${time}/${hist}/normalized.bw \
                --binSize 10 \
                --numberOfProcessors 8
            
            echo "  Created merged: ${bigwigPath}/${org}/${time}/${hist}/normalized.bw"
            
            # Clean up temporary merged BAM
            rm -f ${merged_bam} ${merged_bam}.bai
        fi
    done
done

echo ""
echo "================================================================"
echo "BigWig Normalization Complete!"
echo "================================================================"
echo ""
echo "Output files:"
echo "  Individual replicates: ${bigwigPath}/${org}/{time}/{hist}/normalized_rep{1,2}.bw"
echo "  Merged replicates: ${bigwigPath}/${org}/{time}/{hist}/normalized.bw"
echo ""
echo "Files are RPGC normalized with effective genome size: ${gsize}"
echo "These files are ready for heatmap generation"
echo ""
echo "Total files created:"
echo "  - 24 individual replicate BigWigs (4 conditions × 3 targets × 2 replicates)"
echo "  - 12 merged BigWigs (4 conditions × 3 targets)"
echo "  - IgG samples excluded as requested"