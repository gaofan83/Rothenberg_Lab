mkdir -p juicebox_Bcl11bKO1_DN2
/usr/local/bin//HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i /home/fgao/Data_hiC/Rothenberg_Lab/hic_out_Bcl11bKO1_DN2/hic_results/data/Bcl11bKO1_DN2/Bcl11bKO1_DN2.allValidPairs \
	-g mm10 -j /home/fgao/software/juicer_tools.jar -o /home/fgao/Data_hiC/Rothenberg_Lab/juicebox_Bcl11bKO1_DN2 -r Arima_refrag_mm10.bed

mkdir -p juicebox_Bcl11bKO2_DN2
/usr/local/bin//HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i /home/fgao/Data_hiC/Rothenberg_Lab/hic_out/hic_results/data/Bcl11bKO2_DN2/Bcl11bKO2_DN2.allValidPairs \
        -g mm10 -j /home/fgao/software/juicer_tools.jar -o juicebox_Bcl11bKO2_DN2 -r Arima_refrag_mm10.bed

mkdir -p juicebox_WT1_DN2
/usr/local/bin//HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i /home/fgao/Data_hiC/Rothenberg_Lab/hic_out/hic_results/data/WT1_DN2/WT1_DN2.allValidPairs \
        -g mm10 -j /home/fgao/software/juicer_tools.jar -o juicebox_WT1_DN2 -r Arima_refrag_mm10.bed

mkdir -p juicebox_WT2_DN2
/usr/local/bin//HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i /home/fgao/Data_hiC/Rothenberg_Lab/hic_out/hic_results/data/WT2_DN2/WT2_DN2.allValidPairs \
        -g mm10 -j /home/fgao/software/juicer_tools.jar -o juicebox_WT2_DN2 -r Arima_refrag_mm10.bed
