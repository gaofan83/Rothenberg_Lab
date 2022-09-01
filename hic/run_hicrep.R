library(hicrep)

bin=as.integer(100000)
#mat1 <- cool2matrix(paste("matrix_cool/Bcl11bKO1_DN2_", bin, ".cool", sep=""))
#mat2 <- cool2matrix(paste("matrix_cool/Bcl11bKO2_DN2_", bin, ".cool", sep=""))
#h_value <- htrain(mat1, mat2, resol = bin, lbr = 0, ubr = 5000000, range = 0:10)
h_value <- 3

all.scc <- data.frame(chr=character(), scc=double(), stringsAsFactors=FALSE)
for (i in paste("chr", c(as.character(1:19)), sep="")){
	mat1.chr <- cool2matrix(paste("matrix_cool/Bcl11bKO1_DN2_", bin, ".cool", sep=""), chr = i)
	mat2.chr <- cool2matrix(paste("matrix_cool/Bcl11bKO2_DN2_", bin, ".cool", sep=""), chr = i)
	scc = get.scc(mat1.chr, mat2.chr, bin, h_value, lbr = 0, ubr = 5000000)
	scc_chr <- data.frame(i, scc$scc)
	names(scc_chr) <- c("chr", "scc")
	all.scc<-rbind(all.scc, scc_chr)
}
write.table(all.scc, file="hicrep_Bcl11bKO1_KO2_100kb.txt", quote=F, sep="\t", row.names=F)

all.scc <- data.frame(chr=character(), scc=double(), stringsAsFactors=FALSE)
for (i in paste("chr", c(as.character(1:19)), sep="")){
        mat1.chr <- cool2matrix(paste("matrix_cool/WT1_DN2_", bin, ".cool", sep=""), chr = i)
        mat2.chr <- cool2matrix(paste("matrix_cool/WT2_DN2_", bin, ".cool", sep=""), chr = i)
        scc = get.scc(mat1.chr, mat2.chr, bin, h_value, lbr = 0, ubr = 5000000)
        scc_chr <- data.frame(i, scc$scc)
        names(scc_chr) <- c("chr", "scc")
        all.scc<-rbind(all.scc, scc_chr)
}
write.table(all.scc, file="hicrep_WT1_WT2_100kb.txt", quote=F, sep="\t", row.names=F)

all.scc <- data.frame(chr=character(), scc=double(), stringsAsFactors=FALSE)
for (i in paste("chr", c(as.character(1:19)), sep="")){
        mat1.chr <- cool2matrix(paste("matrix_cool/Bcl11bKO1_DN2_", bin, ".cool", sep=""), chr = i)
        mat2.chr <- cool2matrix(paste("matrix_cool/WT2_DN2_", bin, ".cool", sep=""), chr = i)
        scc = get.scc(mat1.chr, mat2.chr, bin, h_value, lbr = 0, ubr = 5000000)
        scc_chr <- data.frame(i, scc$scc)
        names(scc_chr) <- c("chr", "scc")
        all.scc<-rbind(all.scc, scc_chr)
}
write.table(all.scc, file="hicrep_Bcl11bKO1_WT2_100kb.txt", quote=F, sep="\t", row.names=F)

all.scc <- data.frame(chr=character(), scc=double(), stringsAsFactors=FALSE)
for (i in paste("chr", c(as.character(1:19)), sep="")){
        mat1.chr <- cool2matrix(paste("matrix_cool/WT1_DN2_", bin, ".cool", sep=""), chr = i)
        mat2.chr <- cool2matrix(paste("matrix_cool/Bcl11bKO2_DN2_", bin, ".cool", sep=""), chr = i)
        scc = get.scc(mat1.chr, mat2.chr, bin, h_value, lbr = 0, ubr = 5000000)
        scc_chr <- data.frame(i, scc$scc)
        names(scc_chr) <- c("chr", "scc")
        all.scc<-rbind(all.scc, scc_chr)
}
write.table(all.scc, file="hicrep_WT1_Bcl11bKO2_100kb.txt", quote=F, sep="\t", row.names=F)
