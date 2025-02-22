### Quality Control ###
echo -e "\n\n\n 222# QC Results!!! \n\n\n"
date
QC_DIR="/home/qc"
mkdir -p "$QC_DIR"
cd "$QC_DIR" || exit

fastqc -t 12 -o "$QC_DIR" /home/fastq/*f*.gz
multiqc "$QC_DIR"

echo -e "\n 222# ALL QC Work Done!!! \n"
date

### HISAT2 Index Preparation ###
echo -e "\n\n\n Preparing HISAT2 Index \n\n\n"
date
extract_exons.py /home/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf > /home/reference/genome.exon
extract_splice_sites.py /home/reference/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf > /home/reference/genome.ss
hisat2-build --ss /home/reference/genome.ss --exon /home/reference/genome.exon \
    /home/reference/chm13v2.0.fa /home/reference/T2T/genome_tran

### Alignment with HISAT2 ###
echo -e "\n\n\n 333# Aligning with HISAT2!!! \n\n\n"
date
ALIGN_DIR="/home/T2T/align"
FLAG_DIR="$ALIGN_DIR/flag"
INDEX='/home/reference/T2T/genome_tran '

mkdir -p "$FLAG_DIR"
cd "$ALIGN_DIR" || exit

for file in /home/T2T/fastq/*_1_val_1*gz; do
    id=$(basename "$file" | cut -d'_' -f1)

    echo "333#  ${id} is being processed with HISAT2!"
    hisat2 -t -p 12 -x "$INDEX" \
        -1 "/home/T2T/fastq/${id}*_1_val_1*gz" \
        -2 "/home/T2T/fastq/${id}*_2_val_2*gz" \
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
COUNT_DIR="/home/T2T/counts"
GTF='/home/reference/T2T/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.gtf'
mkdir -p "$COUNT_DIR"
cd "$COUNT_DIR" || exit

featureCounts -a $GTF -F GTF -o repeat_counts.txt -T 4 -M -O -f -g gene_id -t feature -p -B -C -d 10 "$ALIGN_DIR"/*.bam


multiqc "$COUNT_DIR"

echo -e "\n\n\n ALL WORK DONE!!! \n\n\n"
date