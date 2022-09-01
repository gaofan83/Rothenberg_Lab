#!/bin/bash
resolution="5000 10000 20000 50000 100000"
for bin in $resolution; do
	hicConvertFormat -m hic_out_Bcl11bKO1_DN2/hic_results/matrix/Bcl11bKO1_DN2/raw/${bin}/Bcl11bKO1_DN2_${bin}.matrix \
		--bedFileHicpro hic_out_Bcl11bKO1_DN2/hic_results/matrix/Bcl11bKO1_DN2/raw/${bin}/Bcl11bKO1_DN2_${bin}_abs.bed \
		--inputFormat hicpro --outputFormat cool -o Bcl11bKO1_DN2_${bin}.cool
done

resolution="5000 10000 20000 50000 100000"
for bin in $resolution; do
        hicConvertFormat -m hic_out/hic_results/matrix/Bcl11bKO2_DN2/raw/${bin}/Bcl11bKO2_DN2_${bin}.matrix \
                --bedFileHicpro hic_out/hic_results/matrix/Bcl11bKO2_DN2/raw/${bin}/Bcl11bKO2_DN2_${bin}_abs.bed \
                --inputFormat hicpro --outputFormat cool -o Bcl11bKO2_DN2_${bin}.cool
done

resolution="5000 10000 20000 50000 100000"
for bin in $resolution; do
        hicConvertFormat -m hic_out/hic_results/matrix/WT1_DN2/raw/${bin}/WT1_DN2_${bin}.matrix \
                --bedFileHicpro hic_out/hic_results/matrix/WT1_DN2/raw/${bin}/WT1_DN2_${bin}_abs.bed \
                --inputFormat hicpro --outputFormat cool -o WT1_DN2_${bin}.cool
done

resolution="5000 10000 20000 50000 100000"
for bin in $resolution; do
        hicConvertFormat -m hic_out/hic_results/matrix/WT2_DN2/raw/${bin}/WT2_DN2_${bin}.matrix \
                --bedFileHicpro hic_out/hic_results/matrix/WT2_DN2/raw/${bin}/WT2_DN2_${bin}_abs.bed \
                --inputFormat hicpro --outputFormat cool -o WT2_DN2_${bin}.cool
done
