write_RaxML_tree<-function(bootstraps, alnmnt) {
  require("seqinr")
  require("ape")
  require("ips")
  bootstraps->n
  alnmnt<-as.DNAbin(alnmnt)
  raxml(alnmnt, m="GTRGAMMA", N=n, p=5, f="a", x=2, k=FALSE, exec="~/standard-RAxML/raxmlHPC-SSE3", backbone=NULL)
}