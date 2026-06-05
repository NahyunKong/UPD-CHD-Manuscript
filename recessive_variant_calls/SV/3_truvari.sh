LSF_DOCKER_VOLUMES="/storage1/fs1/jin810-2/Active:/storage1/fs1/jin810-2/Active /storage1/fs1/jin810/Active:/storage1/fs1/jin810/Active /storage2/fs1/epigenome/Active/:/storage2/fs1/epigenome/Active/ /home/k.nahyun:/home/k.nahyun" 
bsub -Is -G compute-jin810-t3 -q subscription -sla jin810_t3 -R 'rusage[mem=30GB]' -M 30GB -a 'docker(meredith705/truvari:latest)'  /bin/bash
#!/bin/bash

# Reference genome
ref="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/input_files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"

# Directories
input_dir="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/sv_analysis/SV_calling/result/filtering"   # <-- Update this path
output_dir="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/sv_analysis/SV_calling/result/truvari"  # <-- Update this path
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
    child_hom_vcf="${input_dir}/${sample}_hom_${chr}.vcf.gz"
    mother_het_vcf="${input_dir}/${sample}_mom_het_${chr}.vcf.gz"

    echo "Processing $sample on $chr ..."

    # Truvari: child vs mother
    out_mom="${output_dir}/${sample}_hom_x_mom_het"
    truvari bench -b "$child_hom_vcf" -c "$mother_het_vcf" -f "$ref" \
        --bSample "$sample" --cSample "${sample}-01" -o "$out_mom"

    # Truvari: child vs father
    out_dad="${output_dir}/${sample}_hom_x_father"
    truvari bench -b "$child_hom_vcf" -c "$father_vcf" -f "$ref" \
        --bSample "$sample" --cSample "${sample}-02" -o "$out_dad"

    # Truvari: (hom ∩ mom het) ∩ not in father
    base2="${out_mom}/tp-base.vcf.gz"
    comp3="${out_dad}/fn.vcf.gz"
    out_final="${output_dir}/${sample}_hom_momhet_notfather"
    truvari bench -b "$base2" -c "$comp3" -f "$ref" \
        --bSample "$sample" --cSample "$sample" -o "$out_final"

    echo "Done with $sample"
done

