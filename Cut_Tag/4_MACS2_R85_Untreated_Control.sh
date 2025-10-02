#!/bin/bash
#SBATCH --job-name=MACS2_Separate_Reps_Tutorial
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --output=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%j.out
#SBATCH --error=/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise/logs/%x_%j.err
#SBATCH --mail-user=ayeh.sadr@icr.ac.uk
#SBATCH --mail-type=ALL

# Load modules
module load anaconda/3
source /opt/software/applications/anaconda/3/etc/profile.d/conda.sh
conda activate macs2.2.7.1


# Define paths - following tutorial structure exactly
projPath="/data/scratch/DMP/DUDMP/PAEDCANC/asadr/Cut_Tag_Noise"
bamPath="${projPath}/bam_files_dedup"

# Define experimental variables - following tutorial exactly
org="human"  # Only human samples
genome="hs"
timeControl="Untreated"  # T1 equivalent
TimeCompare=("D7_Cis" "D14_Cis" "D14_C70")  
histList=("H3K4me3" "KDM5A" "KDM5B")

# Sample mapping for each replicate separately
declare -A sample_map
# Replicate 1
sample_map["Untreated_H3K4me3_1"]="R85_0003"
sample_map["D7_Cis_H3K4me3_1"]="R85_0004"
sample_map["D14_Cis_H3K4me3_1"]="R85_0005"
sample_map["D14_C70_H3K4me3_1"]="R85_0006"
sample_map["Untreated_KDM5A_1"]="R85_0015"
sample_map["D7_Cis_KDM5A_1"]="R85_0016"
sample_map["D14_Cis_KDM5A_1"]="R85_0017"
sample_map["D14_C70_KDM5A_1"]="R85_0018"
sample_map["Untreated_KDM5B_1"]="R85_0027"
sample_map["D7_Cis_KDM5B_1"]="R85_0028"
sample_map["D14_Cis_KDM5B_1"]="R85_0029"
sample_map["D14_C70_KDM5B_1"]="R85_0030"

# Replicate 2
sample_map["Untreated_H3K4me3_2"]="R85_0009"
sample_map["D7_Cis_H3K4me3_2"]="R85_0010"
sample_map["D14_Cis_H3K4me3_2"]="R85_0011"
sample_map["D14_C70_H3K4me3_2"]="R85_0012"
sample_map["Untreated_KDM5A_2"]="R85_0021"
sample_map["D7_Cis_KDM5A_2"]="R85_0022"
sample_map["D14_Cis_KDM5A_2"]="R85_0023"
sample_map["D14_C70_KDM5A_2"]="R85_0024"
sample_map["Untreated_KDM5B_2"]="R85_0033"
sample_map["D7_Cis_KDM5B_2"]="R85_0034"
sample_map["D14_Cis_KDM5B_2"]="R85_0035"
sample_map["D14_C70_KDM5B_2"]="R85_0036"

# ======================================
# MAIN PEAK CALLING - FOLLOWING TUTORIAL
# ======================================
echo "Starting MACS2 peak calling following tutorial structure..."
echo "Processing each replicate separately"
echo ""

# Process each replicate separately - following tutorial loop structure
for rep in 1 2; do
    echo "================================================================"
    echo "Processing Replicate ${rep}"
    echo "================================================================"
    
    for time in ${TimeCompare[@]}; do
        for hist in ${histList[@]}; do
            echo "Processing: ${time} ${hist} rep${rep} vs ${timeControl} ${hist} rep${rep}"
            
            # Create output directory - following tutorial path structure
            mkdir -p ${projPath}/peakCalling/MACS2/${org}/${timeControl}_control/${time}/${hist}/rep${rep}
            
            # Get BAM files
            treatment_bam="${bamPath}/${sample_map[${time}_${hist}_${rep}]}_dedup.bam"
            control_bam="${bamPath}/${sample_map[${timeControl}_${hist}_${rep}]}_dedup.bam"
            
            # Check if files exist
            if [ -f "$treatment_bam" ] && [ -f "$control_bam" ]; then
                # Run MACS2 - following tutorial parameters exactly
                macs2 callpeak \
                    -t $treatment_bam \
                    -c $control_bam \
                    -g ${genome} \
                    -f BAMPE \
                    -n macs2_peak_q0.1 \
                    --outdir $projPath/peakCalling/MACS2/${org}/${timeControl}_control/${time}/${hist}/rep${rep}/ \
                    -q 0.1 \
                    --keep-dup all \
                    2>${projPath}/peakCalling/MACS2/${org}/${timeControl}_control/${time}/${hist}/rep${rep}/macs2Peak_summary.txt
            else
                echo "  WARNING: Missing BAM files"
                echo "    Treatment: $treatment_bam"
                echo "    Control: $control_bam"
            fi
        done
    done
done

# ======================================
# RENAME FILES FOR BROWSER - FOLLOWING TUTORIAL
# ======================================
echo ""
echo "Renaming narrowPeak files to .bed format..."

for rep in 1 2; do
    for time in ${TimeCompare[@]}; do
        for hist in ${histList[@]}; do
            peak_file="$projPath/peakCalling/MACS2/${org}/${timeControl}_control/${time}/${hist}/rep${rep}/macs2_peak_q0.1_peaks.narrowPeak"
            if [ -f "$peak_file" ]; then
                mv $peak_file ${peak_file}.bed
                echo "  Renamed: rep${rep} ${time} ${hist}"
            fi
        done
    done
done

# ======================================
# PARAMETER TESTING - FOLLOWING TUTORIAL
# ======================================
echo ""
echo "Testing different parameters for MACS2..."

qValues=(0.1 0.001 0.00001)

for rep in 1 2; do
    echo ""
    echo "Parameter testing for Replicate ${rep}"
    echo "======================================"
    
    for hist in ${histList[@]}; do
        for time in ${TimeCompare[@]}; do
            for q in ${qValues[@]}; do
                echo "Testing: ${time} ${hist} rep${rep} with q=${q}"
                
                # Create directory - note tutorial has typo "mkdir -p outdir", fixing it
                mkdir -p $projPath/peakCalling/MACS2/${org}/parameter/${timeControl}_control/${time}/${hist}/rep${rep}/
                
                # Get BAM files
                treatment_bam="${bamPath}/${sample_map[${time}_${hist}_${rep}]}_dedup.bam"
                control_bam="${bamPath}/${sample_map[${timeControl}_${hist}_${rep}]}_dedup.bam"
                
                if [ -f "$treatment_bam" ] && [ -f "$control_bam" ]; then
                    
                    # normal - following tutorial exactly
                    time macs2 callpeak \
                        -t $treatment_bam \
                        -c $control_bam \
                        -g ${genome} \
                        -f BAMPE \
                        -n macs2_peak_q${q} \
                        --outdir $projPath/peakCalling/MACS2/${org}/parameter/${timeControl}_control/${time}/${hist}/rep${rep}/ \
                        -q ${q} \
                        --keep-dup all
                    
                    # nomodel - following tutorial exactly
                    time macs2 callpeak \
                        -t $treatment_bam \
                        -c $control_bam \
                        -g ${genome} \
                        -f BAMPE \
                        -n macs2_peak_q${q}_nomodel \
                        --outdir $projPath/peakCalling/MACS2/${org}/parameter/${timeControl}_control/${time}/${hist}/rep${rep}/ \
                        -q ${q} \
                        --nomodel \
                        --keep-dup all
                    
                    # broad - following tutorial exactly
                    time macs2 callpeak \
                        -t $treatment_bam \
                        -c $control_bam \
                        -g ${genome} \
                        -f BAMPE \
                        -n macs2_peak_q${q}_broad \
                        --outdir $projPath/peakCalling/MACS2/${org}/parameter/${timeControl}_control/${time}/${hist}/rep${rep}/ \
                        -q ${q} \
                        --broad \
                        --keep-dup all
                    
                    # nomodel and broad - following tutorial exactly
                    time macs2 callpeak \
                        -t $treatment_bam \
                        -c $control_bam \
                        -g ${genome} \
                        -f BAMPE \
                        -n macs2_peak_q${q}_nomodel_broad \
                        --outdir $projPath/peakCalling/MACS2/${org}/parameter/${timeControl}_control/${time}/${hist}/rep${rep}/ \
                        -q ${q} \
                        --broad \
                        --nomodel \
                        --keep-dup all
                fi
            done
        done
    done
done

# ======================================
# SUMMARY REPORT
# ======================================
echo ""
echo "================================================================"
echo "Generating summary report..."
echo "================================================================"

summary_file="${projPath}/peakCalling/MACS2/peak_summary_separate_reps.txt"
echo -e "Replicate\tComparison\tTarget\tQ_value\tMode\tPeak_Count\tSample_IDs" > $summary_file

for rep in 1 2; do
    for time in ${TimeCompare[@]}; do
        for hist in ${histList[@]}; do
            # Get sample IDs for reference
            treatment_id="${sample_map[${time}_${hist}_${rep}]}"
            control_id="${sample_map[${timeControl}_${hist}_${rep}]}"
            
            # Count peaks from main analysis (q=0.1)
            peak_file="${projPath}/peakCalling/MACS2/${org}/${timeControl}_control/${time}/${hist}/rep${rep}/macs2_peak_q0.1_peaks.narrowPeak.bed"
            if [ -f "$peak_file" ]; then
                peak_count=$(wc -l < "$peak_file")
                echo -e "Rep${rep}\t${time}_vs_${timeControl}\t${hist}\t0.1\tdefault\t${peak_count}\t${treatment_id}_vs_${control_id}" >> $summary_file
            fi
            
            # Count peaks from parameter testing
            for q in ${qValues[@]}; do
                for mode in "" "_nomodel" "_broad" "_nomodel_broad"; do
                    if [[ "$mode" == *"broad"* ]]; then
                        ext="broadPeak"
                    else
                        ext="narrowPeak"
                    fi
                    
                    peak_file="${projPath}/peakCalling/MACS2/${org}/parameter/${timeControl}_control/${time}/${hist}/rep${rep}/macs2_peak_q${q}${mode}_peaks.${ext}"
                    if [ -f "$peak_file" ]; then
                        peak_count=$(wc -l < "$peak_file")
                        mode_name=${mode:-"default"}
                        mode_name=${mode_name#_}
                        echo -e "Rep${rep}\t${time}_vs_${timeControl}\t${hist}\t${q}\t${mode_name}\t${peak_count}\t${treatment_id}_vs_${control_id}" >> $summary_file
                    fi
                done
            done
        done
    done
done

echo ""
echo "Peak count summary (q=0.1, default mode):"
echo "=========================================="
grep -E "0.1.*default" $summary_file | column -t

echo ""
echo "================================================================"
echo "MACS2 Peak Calling Complete - Following Tutorial Structure"
echo "================================================================"
echo ""
echo "Directory structure (following tutorial):"
echo "  ${projPath}/peakCalling/MACS2/${org}/"
echo "    ├── ${timeControl}_control/"
echo "    │   ├── D7_Cis/{H3K4me3,KDM5A,KDM5B}/rep{1,2}/"
echo "    │   ├── D14_Cis/{H3K4me3,KDM5A,KDM5B}/rep{1,2}/"
echo "    │   └── D14_C70/{H3K4me3,KDM5A,KDM5B}/rep{1,2}/"
echo "    └── parameter/"
echo "        └── ${timeControl}_control/{time}/{hist}/rep{1,2}/"
echo ""
echo "Each replicate processed separately as per tutorial"
echo "Summary report: ${summary_file}"