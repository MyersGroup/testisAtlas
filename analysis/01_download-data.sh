#!/bin/bash

secret= # put access key here

# If you're reading this post publication, the files will now be on GEO instead

###############
### 10X Hormad1 ###
###############

# use no parent though ????

wget https://htcf.wustl.edu/files/${secret}/MouseHormad1/10X/ -r -nH --cut-dirs=4


###############
### CNP ### (4)
###############

wget https://htcf.wustl.edu/files/${secret}/MouseCNP/010517CNP/run_2131_s_2_1_withindex_sequence.txt_CGTACTAG.fq.gz -O CNP_2017-01-05_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseCNP/010517CNP/run_2131_s_2_3_withindex_sequence.txt_CGTACTAG.fq.gz -O CNP_2017-01-05_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseCNP/121516CNP/run_2131_s_2_1_withindex_sequence.txt_TAAGGCGA.fq.gz -O CNP_2016-12-15_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseCNP/121516CNP/run_2131_s_2_3_withindex_sequence.txt_TAAGGCGA.fq.gz -O CNP_2016-12-15_B.fq.gz


###############
### CUL4A ### (4)
###############

# Cul4a was sequenced twice
#https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0328/run_2192_s_7_1_withindex_sequence.txt_TAAGGCGA.fq.gz
#https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0328/run_2192_s_7_3_withindex_sequence.txt_TAAGGCGA.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0328/run_2192_s_8_1_withindex_sequence.txt_TAAGGCGA.fq.gz -O Cul4a_2017-03-28_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0328/run_2192_s_8_3_withindex_sequence.txt_TAAGGCGA.fq.gz -O Cul4a_2017-03-28_B.fq.gz

#https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0330/run_2192_s_7_1_withindex_sequence.txt_GCTACGCT.fq.gz
#https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0330/run_2192_s_7_3_withindex_sequence.txt_GCTACGCT.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0330/run_2192_s_8_1_withindex_sequence.txt_GCTACGCT.fq.gz -O Cul4a_2017-03-30_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseCul4aKO/0330/run_2192_s_8_3_withindex_sequence.txt_GCTACGCT.fq.gz -O Cul4a_2017-03-30_B.fq.gz


###############
### HORMAD1 ## (6)
###############

wget https://htcf.wustl.edu/files/${secret}/MouseHormad1/Dropseq/run_2348_s_1_1_withindex_sequence.txt_TAAGGCGA.fq.gz -O Hormad1_1_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseHormad1/Dropseq/run_2348_s_1_3_withindex_sequence.txt_TAAGGCGA.fq.gz -O Hormad1_1_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseHormad1/Dropseq/run_2348_s_2_1_withindex_sequence.txt_AGGCAGAA.fq.gz -O Hormad1_2_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseHormad1/Dropseq/run_2348_s_2_3_withindex_sequence.txt_AGGCAGAA.fq.gz -O Hormad1_2_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseHormad1/Dropseq/run_2348_s_2_1_withindex_sequence.txt_TCCTGAGC.fq.gz -O Hormad1_3_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseHormad1/Dropseq/run_2348_s_2_3_withindex_sequence.txt_TCCTGAGC.fq.gz -O Hormad1_3_B.fq.gz


###############
### MLH3 ### (12)
###############

wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/1190Mlh3/s_1_1_withindex_sequence.txt_CGTACTAG.fq.gz -O Mlh3_1190_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/1190Mlh3/s_1_4_withindex_sequence.txt_CGTACTAG.fq.gz -O Mlh3_1190_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/1191Mlh3/s_1_1_withindex_sequence.txt_AGGCAGAA.fq.gz -O Mlh3_1191_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/1191Mlh3/s_1_4_withindex_sequence.txt_AGGCAGAA.fq.gz -O Mlh3_1191_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/1192Mlh3/s_2_1_withindex_sequence.txt_TAAGGCGA.fq.gz -O Mlh3_1192_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/1192Mlh3/s_2_4_withindex_sequence.txt_TAAGGCGA.fq.gz -O Mlh3_1192_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/MLH3KO300treated/s_1_1_withindex_sequence.txt_TCCTGAGC.fq.gz -O Mlh3_300_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/MLH3KO300treated/s_1_4_withindex_sequence.txt_TCCTGAGC.fq.gz -O Mlh3_300_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/MLH3KO800/s_2_1_withindex_sequence.txt_GGACTCCT.fq.gz -O Mlh3_800_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/MLH3KO800/s_2_4_withindex_sequence.txt_GGACTCCT.fq.gz -O Mlh3_800_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/MLH3KO/MLH3-KO-pooled_S1_L001_R1_001.fastq.gz -O Mlh3_pool_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseMLH3KO/MLH3KO/MLH3-KO-pooled_S1_L001_R2_001.fastq.gz -O Mlh3_pool_B.fq.gz

###############
### WILD TYPE ### (28 + 3)
###############

# CHEM (2+3)

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Chem/mj2/1_S1_L001_I1_001.fastq.gz -O mj_2_I.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Chem/mj2/1_S1_L001_R1_001.fastq.gz -O mj_2_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Chem/mj2/1_S1_L001_R2_001.fastq.gz -O mj_2_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Chem/mj3/s_2_1_withindex_sequence.txt_TAAGGCGA.fq.gz -O mj_3_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Chem/mj3/s_2_4_withindex_sequence.txt_TAAGGCGA.fq.gz -O mj_3_B.fq.gz

# FACS (8)

wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPCI/SPCI_S2_L001_R1_001.fastq.gz -O SPCI_1_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPCI/SPCI_S2_L001_R2_001.fastq.gz -O SPCI_1_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPCII/SPCII_S3_L001_R1_001.fastq.gz -O SPCII_1_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPCII/SPCII_S3_L001_R2_001.fastq.gz -O SPCII_1_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPD/SPD_S1_L001_R1_001.fastq.gz -O SPD_1_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPD/SPD_S1_L001_R2_001.fastq.gz -O SPD_1_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPG/run_2131_s_2_1_withindex_sequence.txt_AGGCAGAA.fq.gz -O SPG_1_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/FACS/SPG/run_2131_s_2_3_withindex_sequence.txt_AGGCAGAA.fq.gz -O SPG_1_B.fq.gz

# MEDI (18)

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/3330/28117_S1_L001_R1_001.fastq.gz -O WT_28117_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/3330/28117_S1_L001_R2_001.fastq.gz -O WT_28117_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/110816WT/run_2095_s_2_1_withindex_sequence.txt_AGGCAGAA.fq.gz -O WT_2016-11-08_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/110816WT/run_2095_s_2_3_withindex_sequence.txt_AGGCAGAA.fq.gz -O WT_2016-11-08_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111016WT/run_2096_s_2_1_withindex_sequence.txt_CTCTCTAC.fq.gz -O WT_2016-11-10_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111016WT/run_2096_s_2_3_withindex_sequence.txt_CTCTCTAC.fq.gz -O WT_2016-11-10_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111116WT/run_2096_s_1_1_withindex_sequence.txt_GGACTCCT.fq.gz -O WT_2016-11-11_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111116WT/run_2096_s_1_3_withindex_sequence.txt_GGACTCCT.fq.gz -O WT_2016-11-11_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111516WT/run_2096_s_2_1_withindex_sequence.txt_CAGAGAGG.fq.gz -O WT_2016-11-15_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111516WT/run_2096_s_2_3_withindex_sequence.txt_CAGAGAGG.fq.gz -O WT_2016-11-15_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111616WT/run_2096_s_1_1_withindex_sequence.txt_TAGGCATG.fq.gz -O WT_2016-11-16_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111616WT/run_2096_s_1_3_withindex_sequence.txt_TAGGCATG.fq.gz -O WT_2016-11-16_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111716WT/run_2095_s_1_1_withindex_sequence.txt_CGTACTAG.fq.gz -O WT_2016-11-17_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/111716WT/run_2095_s_1_3_withindex_sequence.txt_CGTACTAG.fq.gz -O WT_2016-11-17_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/112116WT/run_2095_s_1_1_withindex_sequence.txt_TAAGGCGA.fq.gz -O WT_2016-11-21_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/112116WT/run_2095_s_1_3_withindex_sequence.txt_TAAGGCGA.fq.gz -O WT_2016-11-21_B.fq.gz

wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/112216WT/run_2095_s_2_1_withindex_sequence.txt_TCCTGAGC.fq.gz -O WT_2016-11-22_A.fq.gz
wget https://htcf.wustl.edu/files/${secret}/MouseWT/Medi/112216WT/run_2095_s_2_3_withindex_sequence.txt_TCCTGAGC.fq.gz -O WT_2016-11-22_B.fq.gz


# find . -type f -name "*.out" -exec md5sum "{}" + > MD5.chk

# find . -type f -name "*.out" -exec md5sum {} \;

# -exec md5sum {} \;>> /checksums_backup.md5