/usr/local/bin//HiC-Pro_3.1.0/bin/utils/digest_genome.py -r ^GATC G^ANTC -o Arima_refrag_mm10.bed \
    /home/fgao/reference_genome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa
sudo cp /home/fgao/Data_hiC/Rothenberg_Lab/Arima_refrag_mm10.bed /usr/local/bin//HiC-Pro_3.1.0/annotation/
/usr/local/bin//HiC-Pro_3.1.0/bin/HiC-Pro -i /home/fgao/Data_hiC/Rothenberg_Lab/hic_in -o /home/fgao/Data_hiC/Rothenberg_Lab/hic_out -c config-hicpro.txt
