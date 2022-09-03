#!/bin/bash

hicConvertFormat -m matrix_cool/Bcl11bKO1_DN2_5000.cool -o matrix_cool/Bcl11bKO1_DN2_ginteractions \
                --inputFormat cool --outputFormat ginteractions

hicConvertFormat -m matrix_cool/Bcl11bKO2_DN2_5000.cool -o matrix_cool/Bcl11bKO2_DN2_ginteractions \
                --inputFormat cool --outputFormat ginteractions

hicConvertFormat -m matrix_cool/WT1_DN2_5000.cool -o matrix_cool/WT1_DN2_ginteractions \
                --inputFormat cool --outputFormat ginteractions

hicConvertFormat -m matrix_cool/WT2_DN2_5000.cool -o matrix_cool/WT2_DN2_ginteractions \
                --inputFormat cool --outputFormat ginteractions

awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' matrix_cool/Bcl11bKO1_DN2_ginteractions.tsv > matrix_cool/Bcl11bKO1_DN2_ginteractions.tsv.short
sort -k2,2d -k6,6d matrix_cool/Bcl11bKO1_DN2_ginteractions.tsv.short > matrix_cool/Bcl11bKO1_DN2_ginteractions.tsv.short.sorted
~/software/juicer-1.6/scripts/common/juicer_tools pre -r 10000,20000,50000,100000,250000,500000,1000000 \
 matrix_cool/Bcl11bKO1_DN2_ginteractions.tsv.short.sorted Bcl11bKO1_DN2.hic chrom_mm10.sizes

awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' matrix_cool/Bcl11bKO2_DN2_ginteractions.tsv > matrix_cool/Bcl11bKO2_DN2_ginteractions.tsv.short
sort -k2,2d -k6,6d matrix_cool/Bcl11bKO2_DN2_ginteractions.tsv.short > matrix_cool/Bcl11bKO2_DN2_ginteractions.tsv.short.sorted
~/software/juicer-1.6/scripts/common/juicer_tools pre -r 10000,20000,50000,100000,250000,500000,1000000 \
 matrix_cool/Bcl11bKO2_DN2_ginteractions.tsv.short.sorted Bcl11bKO2_DN2.hic chrom_mm10.sizes
 
awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' matrix_cool/WT1_DN2_ginteractions.tsv > matrix_cool/WT1_DN2_ginteractions.tsv.short
sort -k2,2d -k6,6d matrix_cool/WT1_DN2_ginteractions.tsv.short > matrix_cool/WT1_DN2_ginteractions.tsv.short.sorted
~/software/juicer-1.6/scripts/common/juicer_tools pre -r 10000,20000,50000,100000,250000,500000,1000000 \
 matrix_cool/WT1_DN2_ginteractions.tsv.short.sorted WT1_DN2.hic chrom_mm10.sizes

awk -F "\t" '{print 0, $1, $2, 0, 0, $4, $5, 1, $7}' matrix_cool/WT2_DN2_ginteractions.tsv > matrix_cool/WT2_DN2_ginteractions.tsv.short
sort -k2,2d -k6,6d matrix_cool/WT2_DN2_ginteractions.tsv.short > matrix_cool/WT2_DN2_ginteractions.tsv.short.sorted
~/software/juicer-1.6/scripts/common/juicer_tools pre -r 10000,20000,50000,100000,250000,500000,1000000 \
 matrix_cool/WT2_DN2_ginteractions.tsv.short.sorted WT2_DN2.hic chrom_mm10.sizes

#track type=hic name="Bcl11bKO1" bigDataUrl=http://bioinformatics.caltech.edu/events/Ellen/Bcl11bKO1_DN2.hic
