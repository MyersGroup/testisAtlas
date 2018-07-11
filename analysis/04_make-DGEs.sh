#!/bin/bash

###############
### CNP ### (2)
###############

sh 03_generate_DGE.sh CNP_2017-01-05 CNP_2017-01-05_A.fq.gz CNP_2017-01-05_B.fq.gz 1000 >> CNP_2017-01-05_log.txt 2>&1

sh 03_generate_DGE.sh CNP_2016-12-15 CNP_2016-12-15_A.fq.gz CNP_2016-12-15_B.fq.gz 2000 >> CNP_2016-12-15_log.txt 2>&1


###############
### CUL4A ### (2)
###############

sh 03_generate_DGE.sh Cul4a_2017-03-28 run_2192_s_8_1_withindex_sequence.txt_TAAGGCGA.fq.gz run_2192_s_8_3_withindex_sequence.txt_TAAGGCGA.fq.gz >> Cul4a_2017-03-28_log.txt 2>&1

sh 03_generate_DGE.sh Cul4a_2017-03-30 Cul4a_2017-03-30_A.fq.gz Cul4a_2017-03-30_B.fq.gz 3000 >> Cul4a_2017-03-30_log.txt 2>&1

###############
### HORMAD1 ## (3)
###############

sh 03_generate_DGE.sh Hormad1_1_A Hormad1_1_A.fq.gz Hormad1_1_B.fq.gz 3000 >> Hormad1_1_log.txt 2>&1

sh 03_generate_DGE.sh Hormad1_2_A Hormad1_2_A.fq.gz Hormad1_2_B.fq.gz 3000 >> Hormad1_2_log.txt 2>&1

sh 03_generate_DGE.sh Hormad1_3_A Hormad1_3_A.fq.gz Hormad1_3_B.fq.gz 3000 >> Hormad1_3_log.txt 2>&1

###############
### MLH3 ###
###############

sh 03_generate_DGE.sh Mlh3_1190 Mlh3_1190_A.fq.gz Mlh3_1190_B.fq.gz 1500 >> Mlh3_1190_log.txt 2>&1

sh 03_generate_DGE.sh Mlh3_1191 Mlh3_1191_A.fq.gz Mlh3_1191_B.fq.gz 3000 >> Mlh3_1191_log.txt 2>&1

sh 03_generate_DGE.sh Mlh3_1192 Mlh3_1192_A.fq.gz Mlh3_1192_B.fq.gz 3500 >> Mlh3_1192_log.txt 2>&1

sh 03_generate_DGE.sh Mlh3_300 Mlh3_300_A.fq.gz Mlh3_300_B.fq.gz 300 >> Mlh3_300_log.txt 2>&1

sh 03_generate_DGE.sh Mlh3_800 Mlh3_800_A.fq.gz Mlh3_800_B.fq.gz 800 >> Mlh3_800_log.txt 2>&1


###############
### WILD TYPE ###
###############

# CHEM (1)

sh 03_generate_DGE.sh mj_3 mj_3_A.fq.gz mj_3_B.fq.gz 500 >> mj_3_log.txt 2>&1

# FACS (4)

sh 03_generate_DGE.sh SPCI_1 SPCI_1_A.fq.gz SPCI_1_B.fq.gz 200 >> SPCI_1_log.txt 2>&1

sh 03_generate_DGE.sh SPCII_1 SPCII_1_A.fq.gz SPCII_1_B.fq.gz 400 >> SPCII_1_log.txt 2>&1

sh 03_generate_DGE.sh SPD_1 SPD_1_A.fq.gz SPD_1_B.fq.gz 300 >> SPD_1_log.txt 2>&1

sh 03_generate_DGE.sh SPG_1 SPG_1_A.fq.gz SPG_1_B.fq.gz 300 >> SPG_1_log.txt 2>&1

# MEDI (9)

sh 03_generate_DGE.sh WT_28117 WT_28117_A.fq.gz WT_28117_B.fq.gz 1500 >> WT_28117_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-08 WT_2016-11-08_A.fq.gz WT_2016-11-08_B.fq.gz 3000 >> WT_2016-11-08_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-10 WT_2016-11-10_A.fq.gz WT_2016-11-10_B.fq.gz 3000 >> WT_2016-11-10_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-11 WT_2016-11-11_A.fq.gz WT_2016-11-11_B.fq.gz 3000 >> WT_2016-11-11_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-15 WT_2016-11-15_A.fq.gz WT_2016-11-15_B.fq.gz 3000 >> WT_2016-11-15_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-16 WT_2016-11-16_A.fq.gz WT_2016-11-16_B.fq.gz 3000 >> WT_2016-11-16_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-17 WT_2016-11-17_A.fq.gz WT_2016-11-17_B.fq.gz 3000 >> WT_2016-11-17_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-21 WT_2016-11-21_A.fq.gz WT_2016-11-21_B.fq.gz 3000 >> WT_2016-11-21_log.txt 2>&1

sh 03_generate_DGE.sh WT_2016-11-22 WT_2016-11-22_A.fq.gz WT_2016-11-22_B.fq.gz 3000 >> WT_2016-11-22_log.txt 2>&1

