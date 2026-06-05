#!/bin/bash

REFERENCE="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/input_files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
OUTPUT_BASE="/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/sv_analysis/SV_calling/result/sniffles"
THREADS=4

# List of BAM files
bam_files=(
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-03820/m84141_241109_043056_s4.hifi_reads.bc2027.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-03820-01/merged_bam/1-03820-01_merged.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-03820-02/m84063_250211_170620_s1.hifi_reads.bc2024.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-06216-01/m84060_241218_180442_s2.hifi_reads.bc2023.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-06216-02/m84060_241218_180442_s2.hifi_reads.bc2019.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-06216/m84141_241109_043056_s4.hifi_reads.bc2024.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-14785/m84043_241106_010215_s4.hifi_reads.bc2015.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-14785-02/merged_bam/1-14785-02_merged.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-14785-01/merged_bam/1-14785-01_merged.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-02917/m84043_241106_010215_s4.hifi_reads.bc2023.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-02917-02/merged_bam/1-02917-02_merged.bam"
  "/storage1/fs1/jin810/Active/testing/Nahyun/UPD_PROJECT/20241118_long_read/data/1-02917-01/merged_bam/1-02917-01_merged.bam"
)

# Loop over each BAM file and run Sniffles
for bam_file in "${bam_files[@]}"
do
  # Extract sample name as the first directory containing a digit-based ID
  sample_name=$(echo "$bam_file" | grep -oE '/1-[^/]+' | head -n 1 | cut -d'/' -f2)
  output_dir="$OUTPUT_BASE"
  mkdir -p "$output_dir"

  output_vcf="$output_dir/${sample_name}_sniffles.vcf"

  echo "Running Sniffles for sample: $sample_name $output_vcf"
  LSF_DOCKER_VOLUMES="/storage1/fs1/jin810/Active:/storage1/fs1/jin810/Active /storage2/fs1/epigenome/Active/:/storage2/fs1/epigenome/Active/ /home/k.nahyun:/home/k.nahyun" 
  bsub  -G compute-jin810 -q general -n 20 -R 'rusage[mem=60GB]' -a 'docker(ztang301/winnowmap:vSniffles2)' \
    /usr/local/bin/sniffles \
      -t $THREADS \
      --reference "$REFERENCE" \
      -i "$bam_file" \
      --sample-id "$sample_name" \
      -v "$output_vcf"
done
