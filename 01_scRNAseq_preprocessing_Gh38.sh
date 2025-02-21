

### Quality Control ###
echo -e "\n\n\n 222# QC Results!!! \n\n\n"
date
QC_DIR="/home/Gh38/qc"
mkdir -p "$QC_DIR"
cd "$QC_DIR" || exit

fastqc -t 12 -o "$QC_DIR" /home/Gh38/fastq/*f*.gz
multiqc "$QC_DIR"

echo -e "\n 222# ALL QC Work Done!!! \n"
date


extract_exons.py /home/Gh38/reference/gencode.v45.chr_patch_hapl_scaff.annotation.gtf > genome.exon
extract_splice_sites.py /home/Gh38/reference/gencode.v45.chr_patch_hapl_scaff.annotation.gtf > genome.ss
hisat2-build /home/Gh38/reference/GRCh38.p13.genome.fa --ss genome.ss --exon genome.exon genome_tran

### Alignment with HISAT2 ###
echo -e "\n\n\n 333# Aligning with HISAT2!!! \n\n\n"
date
ALIGN_DIR="/home/Gh38/align"
FLAG_DIR="$ALIGN_DIR/flag"
INDEX='/home/Gh38/reference/genome_tran'

mkdir -p "$FLAG_DIR"
cd "$ALIGN_DIR" || exit

for file in /home/Gh38/fastq/*_1_val_1*gz; do
    id=$(basename "$file" | cut -d'_' -f1)

    echo "333#  ${id} is being processed with HISAT2!"
    hisat2 -t -p 12 -x "$INDEX" \
        -1 "/home/Gh38/fastq/${id}*_1_val_1*gz" \
        -2 "/home/Gh38/fastq/${id}*_2_val_2*gz" \
        -S "$ALIGN_DIR/${id}.sam"

    echo -e " ${id} converting SAM to sorted BAM and removing SAM file"
    samtools sort -@ 12 -o "$ALIGN_DIR/${id}_sorted.bam" "$ALIGN_DIR/${id}.sam"
    rm "$ALIGN_DIR/${id}.sam"
done

### BAM Indexing & Flagstat ###
echo -e "\n\n\n Running SAMTOOLS Index and Flagstat \n"
date
cd "$FLAG_DIR" || exit

for bam in "$ALIGN_DIR"/*.bam; do
    samtools flagstat -@ 12 "$bam" > "$(basename "$bam" .bam).flagstat"
done
multiqc "$FLAG_DIR"

echo -e "\n\n\n 333# ALL Alignment Work Done!!! \n\n\n"
date

### FeatureCounts ###
echo -e "\n\n\n 444# Running FeatureCounts!!! \n\n\n"
date
COUNT_DIR="/home/Gh38/counts"
GTF='/home/reference/Gh38/gencode.v45.chr_patch_hapl_scaff.annotation.gtf'
mkdir -p "$COUNT_DIR"
cd "$COUNT_DIR" || exit

featureCounts -T 12 -p -a "$GTF" -o counts.txt "$ALIGN_DIR"/*.bam
multiqc "$COUNT_DIR"

echo -e "\n\n\n ALL WORK DONE!!! \n\n\n"
date
