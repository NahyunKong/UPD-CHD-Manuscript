LSF_DOCKER_VOLUMES="/storage1/fs1/jin810-2/Active:/storage1/fs1/jin810-2/Active /storage1/fs1/jin810/Active:/storage1/fs1/jin810/Active /storage2/fs1/epigenome/Active/:/storage2/fs1/epigenome/Active/ /home/k.nahyun:/home/k.nahyun" 
bsub -Is -G compute-jin810 -q general-interactive -R 'rusage[mem=50GB]' -M 50GB -a 'docker(mgibio/bcftools-cwl:1.12)' /bin/bash

#!/bin/bash

# Reference genome
ref="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/input_files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Directories
input_dir="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/sv_analysis/SV_calling/result/sniffles"   # <-- Update this path
output_dir="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/sv_analysis/SV_calling/result/filtering"  # <-- Update this path
mkdir -p "$output_dir"

# Sample list
samples=("1-02917" "1-03820" "1-06216" "1-14785")

for sample in "${samples[@]}"; do
    # Assign chromosome per sample
    if [[ "$sample" == "1-03820" ]]; then
        chr="chr15"
    else
        chr="chr16"
    fi

    proband_vcf="${input_dir}/${sample}_sniffles.vcf.gz"
    mother_vcf="${input_dir}/${sample}-01_sniffles.vcf.gz"
    father_vcf="${input_dir}/${sample}-02_sniffles.vcf.gz"

    # Output filtered files
    child_hom_vcf="${output_dir}/${sample}_hom_${chr}.vcf.gz"
    mother_het_vcf="${output_dir}/${sample}_mom_het_${chr}.vcf.gz"

    echo "Processing $sample on $chr ..."

    # Filter homozygous child SVs
    bcftools view -i 'GT="1/1" | GT="1|1"' -r "$chr" "$proband_vcf" -Oz -o "$child_hom_vcf"

    # Filter heterozygous mother SVs
    bcftools view -i 'GT="0/1"' -r "$chr" "$mother_vcf" -Oz -o "$mother_het_vcf"


    echo "Done with $sample"
done
