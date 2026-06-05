#!/bin/bash

# Core IDs and corresponding bed files
#CORE_IDS=("1-03573") # List for test1
#CORE_IDS=("1-00660" "1-03573" "1-03820" "1-06216" "1-08240" "1-14523" "1-14785" "1-15783") # List for test
CORE_IDS=("1-00180" "1-00660" "1-02917" "1-03573" "1-03820" "1-04975" "1-06216" "1-08240" "1-14523" "1-14785" "1-15162" "1-15783") # List of with all bams
SUFFIXES=("_father" "_mother" "_child") # Suffixes to append to each core ID

# Base directory for scripts
Job="20240617_UPD_joint_call"
BASE_DIR="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/JOBS"

mkdir -p "$BASE_DIR/$Job"

# Loop through all core IDs
for core_id in "${CORE_IDS[@]}"; do
    id="${core_id}_joint_call" # Joint call ID
    # Script file name
    script_file="${BASE_DIR}/${Job}/${id}.bsub"
    out_dir="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}"

    # Create the output directory
    mkdir -p ${out_dir}

    # Create a new job script
    cat << EOF > ${script_file}
#BSUB -n 16
#BSUB -M 1000000000
#BSUB -R 'select[mem>32GB && tmp>32GB] rusage[mem=40GB, tmp=40GB]'
#BSUB -N
#BSUB -u nahyun@wustl.edu
#BSUB -G compute-jin810
#BSUB -J ${BASE_DIR}/${Job}/${id}.bsub
#BSUB -q general
#BSUB -oo ${BASE_DIR}/${Job}/${id}.bsub.log
#BSUB -a 'docker(elle72/glnexus-gatk:vs1)'
export PATH="$PATH":/gatk:/opt/miniconda/envs/gatk/bin:/opt/miniconda/bin:/opt/miniconda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/glnexus:/glnexus/usr/local/bin:/glnexus:/glnexus/usr/local/bin

cd /storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}
ref="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/reference/hs38DH.fa"
HUMANDB="/storage1/fs1/jin810/Active/References/annovar_20191024/humandb"
f_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/haplotype_caller/${core_id}_father/${core_id}_father_BQSR_gatk_hapcaller.vcf"
m_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/haplotype_caller/${core_id}_mother/${core_id}_mother_BQSR_gatk_hapcaller.vcf"
c_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/haplotype_caller/${core_id}_child/${core_id}_child_BQSR_gatk_hapcaller.vcf"
out_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_joint_call.vcf"
out_recal="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_joint_call_recal.vcf"
out_tranches="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_recal.tranches"
out_vqsr="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_joint_call_recal_vqsr.vcf.gz"

VQSR_RESOURCE_OMNI="--resource:omni,known=false,training=true,truth=true,prior=12.0 /storage1/fs1/jin810/Active/References/known_sites/1000G_omni2.5.hg38.vcf.gz"
VQSR_RESOURCE_1000G="--resource:1000g,known=false,training=true,truth=false,prior=10.0 /storage1/fs1/jin810/Active/References/known_sites/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
VQSR_RESOURCE_dbsnp="--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /storage1/fs1/jin810/Active/References/known_sites/resources-broad-hg38-v0-Homo_sapiens_assembly38.dbsnp138.vcf.gz"
VQSR_RESOURCE_hapmap="--resource:hapmap,known=false,training=true,truth=true,prior=15.0 /storage1/fs1/jin810/Active/References/known_sites/hapmap_3.3.hg38.vcf.gz"
VQSR_RESOURCE_mills="--resource:mills,known=false,training=true,truth=true,prior=12.0 /storage1/fs1/jin810/Active/References/known_sites/resources-broad-hg38-v0-Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
out_lm_sp_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_lm_sp.vcf.gz"
out_reID_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_lm_sp_reID.vcf.gz"
out_annot="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_annot"
out_annot_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_annot.hg38_multianno.vcf"
filtered_vcf="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/${core_id}_annot_filtered.vcf.gz"
# Joint call
echo "####JOINT CALL START####"
#rm -r -f /storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20240613_recessive_calling/new_try_hg38_upd_20240928/results/joint_call/${core_id}/GLnexus.DB
#LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so /opt/conda/bin/glnexus_cli --config gatk -m 5 -t 16 \${c_vcf} \${f_vcf} \${m_vcf}  > \${out_bcf}
gatk --java-options "-Xmx16g" CombineGVCFs     -R \${ref}     --variant \${c_vcf}     --variant \${f_vcf}     --variant \${m_vcf}     -O combined.g.vcf
gatk --java-options "-Xmx16g" GenotypeGVCFs     -R \${ref}     -V combined.g.vcf     -O \${out_vcf}

## VQSR
gatk VariantRecalibrator \
-V \${out_vcf} -O \${out_recal} \
--tranches-file \${out_tranches} \
\${VQSR_RESOURCE_OMNI} \${VQSR_RESOURCE_1000G} \${VQSR_RESOURCE_dbsnp} \${VQSR_RESOURCE_hapmap} \${VQSR_RESOURCE_mills} \
-an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -an MQ \
-tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5  -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0  \
--mode BOTH

gatk  ApplyVQSR \
-R \${ref} -V \${out_vcf} \
--recal-file \${out_recal} \
--tranches-file \${out_tranches} \
-O \${out_vqsr} \
--mode BOTH

#left norm
echo "####vcf conversion START####"
bcftools norm -m-both -f \${ref} -O z -o \${out_lm_sp_vcf} \${out_vqsr}

# reset ID
echo "####reest ID START####"
bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%FIRST_ALT' -O z -o \${out_reID_vcf} \${out_vcf}


#annotation
echo "####annovar annotation START####"
/storage1/fs1/jin810/Active/testing/Nahyun/tools/annovar/table_annovar.pl --vcfinput \${out_reID_vcf} \${HUMANDB} -buildver hg38 -out \${out_annot} -remove \
-protocol refGene,genomicSuperDups,clinvar_20210501,avsnp150,esp6500siv2_all,1000g2015aug_all,exac03,exac03nontcga,gnomad_exome,gnomad_genome,gnomad30_genome,cadd16all,dbSNP,kaviar_20150923,dbnsfp41a,bravo_v8,mcap14,revel,dbscsnv11,fathmmxf_coding,fathmmxf_noncoding,ccREs,SEdb_SE,SEdb_TE \
-operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring .

# Filter variants based on read depth (DP), and genotype quality (GQ)
echo "####filter variants ####"
bcftools view -i 'FORMAT/DP[0]>=8 && FORMAT/GQ[0]>=20' -O z \${out_annot_vcf} -o \${filtered_vcf}

# create index
echo "####creat index START####"
bcftools tabix -p vcf \${filtered_vcf}


EOF
    # Provide execute permissions
    chmod +x "$script_file"

    echo "Created script: $script_file"
done
