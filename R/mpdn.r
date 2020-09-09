mpdn<-function(comm, pd, abundance.weighted = TRUE, time.output=FALSE)
{
  if(sum(colnames(comm)!=rownames(pd))>0)
  {
    sp.name=intersect(colnames(comm),rownames(pd))
    comm=comm[,match(sp.name,colnames(comm))]
    pd=pd[match(sp.name,rownames(pd)),match(sp.name,rownames(pd))]
  }
  comt=comm
  if(!abundance.weighted)
  {
    comt[comt>0]=1
    num=rowSums(comt)
  }
  comt=comt/rowSums(comt)
  comt=as.matrix(comt)
  pd=as.matrix(pd)
  comd=comt %*% pd
  res=comd * comt
  res=rowSums(res)
  if(!abundance.weighted)
  {
    res=res*(num/(num-1))
  }
  res
}