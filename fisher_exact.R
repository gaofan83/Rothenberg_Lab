calc_overlap<-function(set1, set2, total)
  {n_tot<-total
   n_set1<-length(set1)
   n_set2<-length(set2)
   n_pp<-length(intersect(set1, set2))
   n_nn<-n_tot-(n_set1+n_set2-n_pp)
   n_pn<-n_set1-n_pp
   n_np<-n_set2-n_pp
   matri1=matrix(c(n_nn, n_pn, n_np, n_pp),nrow=2)
   #stat1=fisher.test(matri1,alternative="two.sided")
   stat1=fisher.test(matri1,alternative="greater")
   result1=stat1$p.value
   result2=n_pp
   result3=paste(intersect(set1, set2), collapse=",")
   result=list(result1, result2, result3)
   return(result)
}
