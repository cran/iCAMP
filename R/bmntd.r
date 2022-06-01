bmntd<-function(comm, pd, abundance.weighted = TRUE,
                exclude.conspecifics = FALSE,time.output=FALSE,
                unit.sum=NULL,spname.check=TRUE,silent=TRUE)
{
  if(spname.check)
  {
    spc=iCAMP::match.name(cn.list = list(comm=comm),both.list = list(pd=pd),silent = silent)
    comm=spc$comm
    pd=spc$pd
  }
  #if(sum(colnames(comm)!=colnames(pd))>0)
  #{
  #  sp.name=colnames(pd)
  # comm=comm[,match(sp.name,colnames(comm))]
  #}
  comt=comm
  comt[comt>0]=1
  if(!abundance.weighted){com.10=comt}
  N=nrow(comm)
  time1=Sys.time()
  min.d=comm[1,]
  if(exclude.conspecifics)
  {
    pd=as.matrix(pd)
    gc()
    diag(pd)=NA
    for(i in 1:N)
    {
      id=(comm[i,]==0)
      if(sum(!id)==0){comt[i,]=0}else{
        min.d=apply(pd[!id,,drop=FALSE],2,min,na.rm=TRUE)
        comt[i,]=min.d
      }
    }
  }else{
    for(i in 1:N)
    {
      id=(comm[i,]==0)
      if(sum(!id)==0){comt[i,]=0}else{
        min.d[!id]=0
        min.d[id]=apply(pd[!id,id,drop=FALSE],2,min)
        comt[i,]=min.d
      }
    }
  }
  time2=Sys.time()
  if(abundance.weighted)
  {
    if(is.null(unit.sum))
    {
      comm.p=comm/rowSums(comm)
    }else{
      comm.p=comm/unit.sum
    }
    comm.p[rowSums(comm)==0,]=0
    time3=Sys.time()
    res=(as.matrix(comt)) %*% (t(comm.p))
    time4=Sys.time()
    res=(res+t(res))/2
  }else{
    res=(as.matrix(comt)) %*% (t(com.10))
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