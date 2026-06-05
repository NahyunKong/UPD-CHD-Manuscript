#!/bin/bash
# Sample-Chromosome Mapping
declare -A SAMPLE_CHROMOSOME=(
    ["1-00660"]="chr4"
    ["1-08240"]="chr4"
    ["1-04975"]="chr9"
    ["1-03573"]="chr8"
    ["1-15783"]="chr12"
    ["1-00180"]="chr14"
    ["1-03820"]="chr15"
    ["1-02917"]="chr16"
    ["1-06216"]="chr16"
    ["1-14523"]="chr16"
    ["1-14785"]="chr16"
    ["1-15162"]="chr19"
)

# Core IDs
#CORE_IDS=("1-00180" "1-00660" "1-03573" "1-03820" "1-06216" "1-08240" "1-14523" "1-14785" "1-15162" "1-15783")
CORE_IDS=("1-02917" "1-04975" "1-00180")
SUFFIXES=("_father" "_mother" "_child") # Suffixes to append to each core ID

# Base directory for scripts
#Job="20240928_UPD_hg38_recessive_calling"
Job="20240928_UPD_wgs_hg38_recessive_calling"
BASE_DIR="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/JOBS"
PATH=$PATH:/opt/miniconda/bin
mkdir -p "$BASE_DIR/$Job"
# Loop through all core IDs
for core_id in "${CORE_IDS[@]}"; do
    # Get the chromosome for the current core_id
    chromosome=${SAMPLE_CHROMOSOME[$core_id]}
    
    # Loop through all suffixes
    for suffix in "${SUFFIXES[@]}"; do
        id="${core_id}${suffix}" # Combine core ID with suffix
        # Script file name
        script_file="${BASE_DIR}/${Job}/${id}_bqsr_hapcaller.bsub"

        # Create a new job script
        cat << EOF > ${script_file}
#BSUB -n 1
#BSUB -M 4000000
#BSUB -R 'select[mem>4000 && tmp>2] rusage[mem=4000, tmp=2] span[hosts=1]'
#BSUB -N
#BSUB -u nahyun@wustl.edu
#BSUB -G compute-jin810
#BSUB -J ${BASE_DIR}/${Job}/${id}_bqsr_hapcaller.sh
#BSUB -q general
#BSUB -oo ${BASE_DIR}/${Job}/${id}_bqsr_hapcaller.sh.log
#BSUB -a 'docker(broadinstitute/gatk:4.4.0.0)'
cd /storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/
ref="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/reference/hs38DH.fa"
bam="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/dna_bam/wgs/${id}.realign.bam"
known_site="/storage1/fs1/jin810/Active/References/GRCh38/KNOWN_SITES/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
out_bam="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/haplotype_caller/wgs/${id}/${id}_realign_BQSR.bam"
out_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/haplotype_caller/wgs/${id}/${id}_BQSR_gatk_hapcaller.vcf"
table="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/haplotype_caller/wgs/${id}/${id}.recal_data.table"
hap_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/haplotype_caller/wgs/${id}/${id}_BQSR_gatk_hapcaller.g.vcf"

mkdir -p ./results/haplotype_caller/wgs/${id}

# Base recalibration
/gatk/gatk BaseRecalibrator \
   -I \${bam} \
   -R \${ref} \
   --known-sites \${known_site} \
   -O \${table}\
   -L ${chromosome}

# BQSR
/gatk/gatk ApplyBQSR --java-options -Xmx30g -R \${ref} \
-I \${bam} --bqsr-recal-file \${table} -O \${out_bam}\
   -L ${chromosome}

# Run Haplotype Caller
/gatk/gatk HaplotypeCaller --java-options -Xmx30g --input \${out_bam} --output \${hap_vcf} \
--reference \${ref} \
--native-pair-hmm-threads 16 \
-ERC BP_RESOLUTION \
--output-mode EMIT_ALL_ACTIVE_SITES \
   -L ${chromosome}
   
EOF
        # Provide execute permissions
        chmod +x "$script_file"

        echo "Created script: $script_file"
    done
done
