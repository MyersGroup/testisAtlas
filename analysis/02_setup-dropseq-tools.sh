# Download software
cd software
wget https://github.com/broadinstitute/picard/releases/download/2.12.0/picard.jar

curl -o dropseq_tools http://mccarrolllab.com/download/922/
unzip dropseq_tools

wget https://github.com/alexdobin/STAR/releases/download/2.5.0b/Linux_x86_64_static_gcc.tgz
tar -zxvf Linux_x86_64_static_gcc.tgz
rm Linux_x86_64_static_gcc.tgz

cd ../


# get up to date metadata
cd metadata

wget ftp://ftp.ensembl.org/pub/release-90/gtf/mus_musculus/Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf.gz
sum Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf.gz
#33985 27861
gunzip Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf.gz

wget ftp://ftp.ensembl.org/pub/release-90/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz
sum Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz
#37031 862254
gunzip Mus_musculus.GRCm38.dna_sm.toplevel.fa.gz

java -jar software/picard.jar CreateSequenceDictionary \
      R=metadata/Mus_musculus.GRCm38.dna_sm.toplevel.fa \
      O=metadata/Mus_musculus.GRCm38.dna_sm.toplevel.dict

sh software/Drop-seq_tools-1.12/ConvertToRefFlat \
	ANNOTATIONS_FILE=metadata/Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf \
	SEQUENCE_DICTIONARY=metadata/Mus_musculus.GRCm38.dna_sm.toplevel.dict \
	OUTPUT=metadata/Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.refFlat

cd ../

# Dropseq Metadata - ensembl release 75?
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63472/suppl/GSE63472_mm10_reference_metadata.tar.gz
#tar -zxvf GSE63472_mm10_reference_metadata.tar.gz
#rm GSE63472_mm10_reference_metadata.tar.gz

# Create STAR genome
software/Linux_x86_64_static_gcc5.3.0/STAR \
	--runMode genomeGenerate \
	--genomeDir metadata \
	--genomeFastaFiles metadata/Mus_musculus.GRCm38.dna_sm.toplevel.fa \
	--runThreadN 5 \
	--sjdbGTFfile metadata/Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf \
	--sjdbOverhang 74

# rename files for use with DropSeq tools
mv metadata/Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.refFlat metadata/Mus_musculus.GRCm38.dna_sm.toplevel.refFlat
mv Mus_musculus.GRCm38.dna_sm.toplevel.fa Mus_musculus.GRCm38.dna_sm.toplevel.fasta