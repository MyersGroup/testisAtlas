#!/bin/bash

SAMPLE=$1
F1=$2
F2=$3
NUM_CELLS=$4

echo $SAMPLE
echo $F1
echo $F2
echo $NUM_CELLS

mkdir results/$SAMPLE


java -Xmx4g -jar software/picard.jar FastqToSam \
      FASTQ=raw_data/$F1 \
	  FASTQ2=raw_data/$F2 \
      OUTPUT=results/$SAMPLE/sam_file.bam \
      SAMPLE_NAME=$SAMPLE


sh software/Drop-seq_tools-1.12/Drop-seq_alignment.sh \
	-g metadata \
	-d software/Drop-seq_tools-1.12 \
	-o results/$SAMPLE \
	-t tmp \
	-r metadata/Mus_musculus.GRCm38.dna_sm.toplevel.fasta \
	-s software/Linux_x86_64_static_gcc5.3.0/STAR \
	-p \
	results/$SAMPLE/sam_file.bam


software/Drop-seq_tools-1.12/DetectBeadSynthesisErrors \
	INPUT=results/$SAMPLE/star_gene_exon_tagged.bam \
	OUTPUT=results/$SAMPLE/star_gene_exon_tagged_corrected.bam \
	OUTPUT_STATS=results/$SAMPLE/my.synthesis_stats.txt \
	SUMMARY=results/$SAMPLE/my.synthesis_stats.summary.txt \
	NUM_BARCODES=$(expr $NUM_CELLS \* 2) \
	PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC


software/Drop-seq_tools-1.12/BAMTagHistogram \
	I=results/$SAMPLE/star_gene_exon_tagged_corrected.bam \
	O=results/$SAMPLE/out_cell_readcounts.txt.gz \
	TAG=XC

#R
#a=read.table("results/out_cell_readcounts.txt.gz", header=F, stringsAsFactors=F)
#x=cumsum(a$V1)
#x=x/max(x)
#
#pdf()
#str(a)
#plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", xlim=c(1,100000))
#dev.off()

# generate DGE
software/Drop-seq_tools-1.12/DigitalExpression \
	I=results/$SAMPLE/star_gene_exon_tagged_corrected.bam \
	O=results/$SAMPLE/out_gene_exon_tagged.dge.txt.gz \
	SUMMARY=results/$SAMPLE/out_gene_exon_tagged.dge.summary.txt \
	NUM_CORE_BARCODES=$NUM_CELLS

# move star logs
mv tmp/star.* results/$SAMPLE/

cd results/$SAMPLE

# Add sample name to output files
for filename in *; do mv "$filename" "${SAMPLE}_$filename"; done;

