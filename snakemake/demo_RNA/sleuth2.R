library(sleuth)
s2c<-read.table("sample_sheet2.txt",header=T,stringsAsFactors=FALSE)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl",
  host = "jul2019.archive.ensembl.org")
ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
  "ensembl_gene_id", "external_gene_name", "description",
  "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))

so1 <- sleuth_prep(s2c, target_mapping=ttg, aggregation_column = 'ens_gene', extra_bootstrap_summary = TRUE)

so1 <- sleuth_fit(so1, ~1, 'reduced')
so1 <- sleuth_fit(so1, ~group, 'full')
so1 <- sleuth_lrt(so1, 'reduced', 'full')

sleuth_table_gene <- sleuth_results(so1, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)
dim(sleuth_table_gene)

sleuth_table_tx <- sleuth_results(so1, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
sleuth_table_tx <- dplyr::filter(sleuth_table_tx, qval <= 0.05)
dim(sleuth_table_tx)

write.table(sleuth_table_gene, file="sleuth_table_gene_deg_group.txt", quote=F,row.names=F,sep="\t")
write.table(sleuth_table_tx, file="sleuth_table_tx_deg_group.txt", quote=F,row.names=F,sep="\t")

