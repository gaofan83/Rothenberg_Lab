library(multiHiCcompare)
data_dir<-"/home/fgao/Data_hiC/Rothenberg_Lab/"

mat_KO1 <- read.table(paste(data_dir, "hic_out_Bcl11bKO1_DN2/hic_results/matrix/Bcl11bKO1_DN2/raw/100000/Bcl11bKO1_DN2_100000.matrix", sep=""))
bed_KO1 <- read.table(paste(data_dir, "hic_out_Bcl11bKO1_DN2/hic_results/matrix/Bcl11bKO1_DN2/raw/100000/Bcl11bKO1_DN2_100000_abs.bed", sep=""))
dat_KO1 <- HiCcompare::hicpro2bedpe(mat_KO1, bed_KO1)

mat_KO2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/Bcl11bKO2_DN2/raw/100000/Bcl11bKO2_DN2_100000.matrix", sep=""))
bed_KO2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/Bcl11bKO2_DN2/raw/100000/Bcl11bKO2_DN2_100000_abs.bed", sep=""))
dat_KO2 <- HiCcompare::hicpro2bedpe(mat_KO2, bed_KO2)

mat_WT1 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT1_DN2/raw/100000/WT1_DN2_100000.matrix", sep=""))
bed_WT1 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT1_DN2/raw/100000/WT1_DN2_100000_abs.bed", sep=""))
dat_WT1 <- HiCcompare::hicpro2bedpe(mat_WT1, bed_WT1)

mat_WT2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT2_DN2/raw/100000/WT2_DN2_100000.matrix", sep=""))
bed_WT2 <- read.table(paste(data_dir, "hic_out/hic_results/matrix/WT2_DN2/raw/100000/WT2_DN2_100000_abs.bed", sep=""))
dat_WT2 <- HiCcompare::hicpro2bedpe(mat_WT2, bed_WT2)

# Install TopDom
# remotes::install_github("HenrikBengtsson/TopDom", ref="master")
library(TopDom)
library(HiCcompare)
library(readr)
library(GenomicRanges)
for (i in 1:(length(dat_KO1$cis)-3))
{
        dat_KO1_chr<-data.frame(dat_KO1$cis[i])
        dat_KO2_chr<-data.frame(dat_KO2$cis[i])
        dat_WT1_chr<-data.frame(dat_WT1$cis[i])
        dat_WT2_chr<-data.frame(dat_WT2$cis[i])
        
        mat_KO1_chr<-data.frame(region1=dat_KO1_chr[,2], region2=dat_KO1_chr[,5], IF=dat_KO1_chr[,7], stringsAsFactors = FALSE)
        mat_KO1_chr_nn <- sparse2full(mat_KO1_chr)
        bed_KO1_chr<-data.frame(chr=rep(dat_KO1_chr[1,1], nrow(mat_KO1_chr_nn)), 
                                    start=as.numeric(rownames(mat_KO1_chr_nn)), end=as.numeric(rownames(mat_KO1_chr_nn))+100000)
        mat_KO1_chr_all <- cbind(bed_KO1_chr, mat_KO1_chr_nn)
        
        mat_KO2_chr<-data.frame(region1=dat_KO2_chr[,2], region2=dat_KO2_chr[,5], IF=dat_KO2_chr[,7], stringsAsFactors = FALSE)
        mat_KO2_chr_nn <- sparse2full(mat_KO2_chr)
        bed_KO2_chr<-data.frame(chr=rep(dat_KO2_chr[1,1], nrow(mat_KO2_chr_nn)), 
                                    start=as.numeric(rownames(mat_KO2_chr_nn)), end=as.numeric(rownames(mat_KO2_chr_nn))+100000)
        mat_KO2_chr_all <- cbind(bed_KO2_chr, mat_KO2_chr_nn)
        
        mat_WT1_chr<-data.frame(region1=dat_WT1_chr[,2], region2=dat_WT1_chr[,5], IF=dat_WT1_chr[,7], stringsAsFactors = FALSE)
        mat_WT1_chr_nn <- sparse2full(mat_WT1_chr)
        bed_WT1_chr<-data.frame(chr=rep(dat_WT1_chr[1,1], nrow(mat_WT1_chr_nn)), 
                                    start=as.numeric(rownames(mat_WT1_chr_nn)), end=as.numeric(rownames(mat_WT1_chr_nn))+100000)
        mat_WT1_chr_all <- cbind(bed_WT1_chr, mat_WT1_chr_nn)
        
        mat_WT2_chr<-data.frame(region1=dat_WT2_chr[,2], region2=dat_WT2_chr[,5], IF=dat_WT2_chr[,7], stringsAsFactors = FALSE)
        mat_WT2_chr_nn <- sparse2full(mat_WT2_chr)
        bed_WT2_chr<-data.frame(chr=rep(dat_WT2_chr[1,1], nrow(mat_WT2_chr_nn)), 
                                    start=as.numeric(rownames(mat_WT2_chr_nn)), end=as.numeric(rownames(mat_WT2_chr_nn))+100000)
        mat_WT2_chr_all <- cbind(bed_WT2_chr, mat_WT2_chr_nn)
        
        write_tsv(mat_KO1_chr_all, paste(data_dir, "KO1_hic_chr", i, ".matrix", sep=""), col_names=FALSE)
        write_tsv(mat_KO2_chr_all, paste(data_dir, "KO2_hic_chr", i, ".matrix", sep=""), col_names=FALSE)
        write_tsv(mat_WT1_chr_all, paste(data_dir, "WT1_hic_chr", i, ".matrix", sep=""), col_names=FALSE)
        write_tsv(mat_WT2_chr_all, paste(data_dir, "WT2_hic_chr", i, ".matrix", sep=""), col_names=FALSE)
        
        TADs_KO1 <- TopDom(paste(data_dir, "KO1_hic_chr", i, ".matrix", sep=""), window.size = 5, outFile = paste(data_dir, "KO1_hic_chr",i, sep=""))
        boundaries_KO1 <- TADs_KO1$bed
        boundaries_KO1 <- boundaries_KO1[boundaries_KO1$name == "boundary",]
        boundaries_KO1 <- makeGRangesFromDataFrame(boundaries_KO1, seqnames.field = 'chrom', start.field = 'chromStart',
                          end.field = 'chromEnd', keep.extra.columns = TRUE)

        TADs_KO2 <- TopDom(paste(data_dir, "KO2_hic_chr", i, ".matrix", sep=""), window.size = 5, outFile = paste(data_dir, "KO2_hic_chr",i, sep=""))
        boundaries_KO2 <- TADs_KO2$bed
        boundaries_KO2 <- boundaries_KO2[boundaries_KO2$name == "boundary",]
        boundaries_KO2 <- makeGRangesFromDataFrame(boundaries_KO2, seqnames.field = 'chrom', start.field = 'chromStart',
                          end.field = 'chromEnd', keep.extra.columns = TRUE)

        TADs_WT1 <- TopDom(paste(data_dir, "WT1_hic_chr", i, ".matrix", sep=""), window.size = 5, outFile = paste(data_dir, "WT1_hic_chr",i, sep=""))
        boundaries_WT1 <- TADs_WT1$bed
        boundaries_WT1 <- boundaries_WT1[boundaries_WT1$name == "boundary",]
        boundaries_WT1 <- makeGRangesFromDataFrame(boundaries_WT1, seqnames.field = 'chrom', start.field = 'chromStart',
                          end.field = 'chromEnd', keep.extra.columns = TRUE)

        TADs_WT2 <- TopDom(paste(data_dir, "WT2_hic_chr", i, ".matrix", sep=""), window.size = 5, outFile = paste(data_dir, "WT2_hic_chr",i, sep=""))
        boundaries_WT2 <- TADs_WT2$bed
        boundaries_WT2 <- boundaries_WT2[boundaries_WT2$name == "boundary",]
        boundaries_WT2 <- makeGRangesFromDataFrame(boundaries_WT2, seqnames.field = 'chrom', start.field = 'chromStart',
                          end.field = 'chromEnd', keep.extra.columns = TRUE)
}

library(tidyverse)
library(pheatmap)
library(svglite)
data_dir<-"/home/fgao/Data_hiC/Rothenberg_Lab/"

for (i in 1:19)
{
  ko1_signal<-read.table(paste(data_dir, "KO1_hic_chr", i, ".binSignal", sep=""), header=T)
  ko2_signal<-read.table(paste(data_dir, "KO2_hic_chr", i, ".binSignal", sep=""), header=T)
  wt1_signal<-read.table(paste(data_dir, "WT1_hic_chr", i, ".binSignal", sep=""), header=T)
  wt2_signal<-read.table(paste(data_dir, "WT2_hic_chr", i, ".binSignal", sep=""), header=T)
  
  ko1_signal$region<-paste(ko1_signal$chr, ":", ko1_signal$from.coord, "-", ko1_signal$to.coord, sep="")
  ko2_signal$region<-paste(ko2_signal$chr, ":", ko2_signal$from.coord, "-", ko2_signal$to.coord, sep="")
  wt1_signal$region<-paste(wt1_signal$chr, ":", wt1_signal$from.coord, "-", wt1_signal$to.coord, sep="")
  wt2_signal$region<-paste(wt2_signal$chr, ":", wt2_signal$from.coord, "-", wt2_signal$to.coord, sep="")
  list_all<-list(ko1_signal,ko2_signal,wt1_signal,wt2_signal)
  list_all_merge <- list_all %>% reduce(inner_join, by='region')
  
  cf_ave<-list_all_merge[,c(6, 14, 21, 28)]
  colnames(cf_ave)<-c("KO1", "KO2", "WT1", "WT2")
  
  cf_ave_cor<-cor(cf_ave, use = "complete.obs")
  pheatmap(cf_ave_cor, main = paste("Correlation of HiC Interact Freq, Chr", i, sep=""), filename=paste(data_dir, "cor_hic_if_chr", i, ".png", sep=""))
  
}
