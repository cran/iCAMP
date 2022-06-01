bmntd.big<-function(comm, pd.desc="pd.desc", pd.spname, pd.wd,
                    spname.check=FALSE, abundance.weighted = TRUE,
                    exclude.conspecifics = FALSE,time.output=FALSE)
{
  requireNamespace("bigmemory")
  
  if(spname.check)
  {
    check.sp=iCAMP::match.name(name.check = pd.spname,cn.list = list(comm=comm))
    comm=check.sp$comm
  }
  pdbig.id=match(colnames(comm),pd.spname)
  comt=comm
  comt[comt>0]=1
  if(!abundance.weighted){com.10=comt}
  N=nrow(comm)
  time1=Sys.time()
  pd=try(bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc))))
  if(inherits(pd,"try-error")){pd=bigmemory::attach.big.matrix(paste0(pd.wd,"/",pd.desc))}
  sp.num=ncol(comm)
  min.d=rep(0,sp.num)
  if(exclude.conspecifics)
  {
    for(i in 1:N)
    {
      id=which(comm[i,]>0)
      pdii=pd[pdbig.id[id],pdbig.id,drop=FALSE]
      pdii[matrix(c(1:nrow(pdii),id),ncol=2)]=NA
      min.d=apply(pdii,2,min,na.rm=TRUE)
      comt[i,]=min.d
    }
  }else{
    for(i in 1:N)
    {
      id=(comm[i,]==0)
      min.d[!id]=0
      min.d[id]=apply(pd[pdbig.id[!id],pdbig.id[id],drop=FALSE],2,min)
      comt[i,]=min.d
    }
  }
  time2=Sys.time()
  if(abundance.weighted)
  {
    comm.p=comm/rowSums(comm)
    time3=Sys.time()
    res=as.matrix(comt) %*% (as.matrix(t(comm.p)))
    time4=Sys.time()
    res=(res+t(res))/2
  }else{
    res=as.matrix(comt) %*% (as.matrix(t(com.10)))
    time3=Sys.time()
    samp.n=rowSums(com.10)
    com.n=matrix(samp.n,nrow = N,ncol = N)
    com.n=com.n+t(com.n)
    time4=Sys.time()
    res=(res+t(res))/com.n
  }
  res=stats::as.dist(res)
  time5=Sys.time()
  if(time.output)
  {
    time=c(time5,time4,time3,time2)-c(time4,time3,time2,time1)
    output=list(result=res,time=time)
  }else{
    output=res
  }
  output
}