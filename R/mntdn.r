mntdn<-function(comm, pd, abundance.weighted = TRUE,
                check.name=TRUE,memory.G=50,time.count=FALSE)
{
  if(time.count) t1=Sys.time()
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memory.G*1024)
    {
      memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
      if(class(memotry)=="try-error"){warning(memotry[1])}
    }
  }
  if(check.name)
  {
    spc=iCAMP::match.name(cn.list=list(comm=comm),both.list=list(pd=pd))
    comm=spc$comm
    pd=spc$pd
  }
  pd=as.matrix(pd)
  N=nrow(comm)
  id=(comm>0)
  diag(pd)=NA
  gc()
  if(abundance.weighted)
  {
    min.d=matrix(0,nrow = N,ncol = ncol(comm))
    for(i in 1:N)
    {
      pdx=pd[id[i,],id[i,],drop=FALSE]
      if(nrow(pdx)>1){min.d[i,id[i,]]=apply(pdx,2,min,na.rm=TRUE)}
    }
    comm.p=comm/rowSums(comm)
    res=min.d * comm.p
    res=rowSums(res)
  }else{
    res=comm[,1]
    for(i in 1:N)
    {
      pdx=pd[id[i,],id[i,],drop=FALSE]
      res[i]=mean(apply(pdx,2,min,na.rm=TRUE))
    }
  }
  if(time.count){t2=format(Sys.time()-t1);message("----Phylogenetic distance calculation costed ",t2," in total----")}
  names(res)=rownames(comm)
  res
}