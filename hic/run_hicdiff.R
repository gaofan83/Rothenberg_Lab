library(multiHiCcompare)
data_dir<-"/home/fgao/Data_hiC/Rothenberg_Lab/"

mat_KO1 <- read.table(paste(data_dir, "hic_out_Bcl11bKO1_DN2/hic_results/matrix/Bcl11bKO1_DN2/raw/10000/Bcl11bKO1_DN2_10000.matrix", sep=""))
bed_KO1 <- read.table(paste(data_dir, "hic_out_Bcl11bKO1_DN2/hic_results/matrix/Bcl11bKO1_DN2/raw/10000/Bcl11bKO1_DN2_10000_abs.bed", sep=""))
dat_KO1 <- HiCcompare::hicpro2bedpe(mat_KO1, bed_KO1)

mat_KO2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/Bcl11bKO2_DN2/raw/10000/Bcl11bKO2_DN2_10000.matrix", sep=""))
bed_KO2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/Bcl11bKO2_DN2/raw/10000/Bcl11bKO2_DN2_10000_abs.bed", sep=""))
dat_KO2 <- HiCcompare::hicpro2bedpe(mat_KO2, bed_KO2)

mat_WT1 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT1_DN2/raw/10000/WT1_DN2_10000.matrix", sep=""))
bed_WT1 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT1_DN2/raw/10000/WT1_DN2_10000_abs.bed", sep=""))
dat_WT1 <- HiCcompare::hicpro2bedpe(mat_WT1, bed_WT1)

mat_WT2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT2_DN2/raw/10000/WT2_DN2_10000.matrix", sep=""))
bed_WT2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT2_DN2/raw/10000/WT2_DN2_10000_abs.bed", sep=""))
dat_WT2 <- HiCcompare::hicpro2bedpe(mat_WT2, bed_WT2)


for (i in 1:(length(dat_KO1$cis)-1))
{
        dat_KO1_chr<-data.frame(chr=data.frame(dat_KO1$cis[i])[,1], region1=(data.frame(dat_KO1$cis[i])[,2] + data.frame(dat_KO1$cis[i])[,3])/2,
                        region2=(data.frame(dat_KO1$cis[i])[,5] + data.frame(dat_KO1$cis[i])[,6])/2, IF=data.frame(dat_KO1$cis[i])[,7])

        dat_KO2_chr<-data.frame(chr=data.frame(dat_KO2$cis[i])[,1], region1=(data.frame(dat_KO2$cis[i])[,2] + data.frame(dat_KO2$cis[i])[,3])/2,
                        region2=(data.frame(dat_KO2$cis[i])[,5] + data.frame(dat_KO2$cis[i])[,6])/2, IF=data.frame(dat_KO2$cis[i])[,7])

        dat_WT1_chr<-data.frame(chr=data.frame(dat_WT1$cis[i])[,1], region1=(data.frame(dat_WT1$cis[i])[,2] + data.frame(dat_WT1$cis[i])[,3])/2,
                        region2=(data.frame(dat_WT1$cis[i])[,5] + data.frame(dat_WT1$cis[i])[,6])/2, IF=data.frame(dat_WT1$cis[i])[,7])

        dat_WT2_chr<-data.frame(chr=data.frame(dat_WT2$cis[i])[,1], region1=(data.frame(dat_WT2$cis[i])[,2] + data.frame(dat_WT2$cis[i])[,3])/2,
                        region2=(data.frame(dat_WT2$cis[i])[,5] + data.frame(dat_WT2$cis[i])[,6])/2, IF=data.frame(dat_WT2$cis[i])[,7])

        hicexp <- make_hicexp(dat_KO1_chr, dat_KO2_chr, dat_WT1_chr, dat_WT2_chr, groups = c(1, 1, 2, 2))

# jointly normalize data
        hicexp <- cyclic_loess(hicexp, parallel =FALSE)

# compare groups
        hicexp <- hic_exactTest(hicexp, parallel = FALSE)

# view manhattan plot of results
#       png(paste(data_dir, "diffhic_manhattan", names(dat_KO1$cis[i]), ".png", sep=""), width=640, height=480)
#       manhattan_hicexp(hicexp)
#       dev.off()

        saveRDS(hicexp,file=paste(data_dir, "hiccompare_", names(dat_KO1$cis[i]), ".rds", sep=""),compress=F)
}


for (i in list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")){
   data<-readRDS(paste("hiccompare_", i, ".rds", sep=""))
   write.table(data@comparison, file=paste("hiccompare_", i, "_stat.txt", sep=""), quote=F, sep="\t", row.names=F)
}

bedgraph<-data.frame(chr=character(), start=integer(), end=integer(), logpaj=double())

for (i in list("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chrX")){
   data<-read.table(paste("hiccompare_", i, "_stat.txt", sep=""), header=T)
   bedgraph_chr<-data.frame(chr=i, start=data$region1, end=data$region2, logpadj=log10(data$p.adj)*(-1))
   bedgraph<-rbind(bedgraph, bedgraph_chr)
}
write.table(bedgraph, file="hiccompare_stat.bedGraph", quote=F, sep="\t", row.names=F, col.names=F)
