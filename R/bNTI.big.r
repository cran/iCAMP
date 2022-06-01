bNTI.big<-function(comm, meta.group=NULL, pd.desc="pd.desc",
                   pd.spname,pd.wd, spname.check=TRUE,
                   nworker=4, memo.size.GB=50, weighted=TRUE,
                   exclude.consp=FALSE,rand=1000,output.dtail=FALSE,
                   RC=FALSE, trace=TRUE)
{
  requireNamespace("parallel")
  requireNamespace("bigmemory")
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memo.size.GB*1024)
    {
      memotry=try(utils::memory.limit(size=memo.size.GB*1024),silent = TRUE)
      if(inherits(memotry,"try-error")){warning(memotry[1])}
    }
  }
  
  if(spname.check)
  {
    aslist<-function(a){if(is.null(a)){NULL}else{out=list(a);names(out)=deparse(substitute(a));out}}
    sampc=iCAMP::match.name(rn.list = c(aslist(comm),aslist(meta.group)))
    comm=sampc$comm
    if(!is.null(meta.group)) meta.group=sampc$meta.group
    check.sp=match.name(name.check = pd.spname,cn.list = list(comm=comm))
    comm=check.sp$comm
  }
  
  samp.name=rownames(comm)
  if(trace){trace.seq=floor(seq(from=1,to=rand,by=50))}else{trace.seq=NULL}
  
  ## randomization function ##
  sp.num=length(pd.spname)
  
  permx<-function(comx,rand)
  {
    nnz=which(colSums(comx)>0)
    permat=lapply(1:rand,
                  function(i)
                  {
                    out=1:ncol(comx)
                    out[nnz]=nnz[sample(1:length(nnz),size=length(nnz))]
                    out
                  })
    permat
  }
  
  if(is.null(meta.group))
  {
    permat=permx(comx=comm,rand=rand)
  }else{
    meta.lev=unique(meta.group[,1])
    permat=list()
    for(j in 1:length(meta.lev))
    {
      sampj=rownames(meta.group)[which(meta.group[,1]==meta.lev[j])]
      permat[[j]]=permx(comx=comm[which(rownames(comm) %in% sampj),,drop=FALSE],rand=rand)
    }
  }
 
  bMNTD.random<-function(i,permat,pd.desc,com,meta.group,pd.spname,
                         weighted,exclude.consp,pd.wd,trace.seq)
  {
    if(!is.null(trace.seq)){if(i %in% trace.seq) message("Now randomizing i=",i,". ",date())}
    
    if(is.null(meta.group))
    {
      comr=com[,permat[[i]],drop=FALSE]
    }else{
      meta.lev=unique(meta.group[,1])
      comr=com
      for(j in 1:length(meta.lev))
      {
        idj=which(rownames(com) %in% (rownames(meta.group)[which(meta.group[,1]==meta.lev[j])]))
        comr[idj,]=com[idj,permat[[j]][[i]]]
      }
    }
    rownames(comr)=rownames(com)
    colnames(comr)=colnames(com)
    bMNTD.rand<-iCAMP::bmntd.big(comm=comr, pd.desc=pd.desc, pd.spname=pd.spname,
                          pd.wd=pd.wd, spname.check=FALSE, abundance.weighted=weighted,
                          exclude.conspecifics=exclude.consp, time.output = FALSE)
    as.matrix(bMNTD.rand)
  }
  
  ####################
  
  if(trace) message("Now calculating observed betaMNTD. Begin at ", date(),". Please wait...")
  bMNTD.obs<-as.matrix(bmntd.big(comm=comm, pd.desc=pd.desc, pd.spname=pd.spname, pd.wd=pd.wd,
                                 spname.check=FALSE, abundance.weighted=weighted,
                                 exclude.conspecifics=exclude.consp,time.output=FALSE))
  gc()
  
  if(nworker==1)
  {
    bMNTD.rand=lapply(1:rand,bMNTD.random,permat=permat,pd.desc=pd.desc,
                      com=comm,meta.group=meta.group,pd.spname=pd.spname,
                      weighted=weighted,exclude.consp=exclude.consp,
                      pd.wd=pd.wd,trace.seq=trace.seq)
  }else{
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
    if(trace) message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
    bMNTD.rand<-parallel::parLapply(c1,1:rand,bMNTD.random,permat=permat,pd.desc=pd.desc,
                                    com=comm,meta.group=meta.group,pd.spname=pd.spname,
                                    weighted=weighted,exclude.consp=exclude.consp,
                                    pd.wd=pd.wd,trace.seq=NULL)
    parallel::stopCluster(c1)
    gc()
    message("Now parallel computation has finished, beta NTI is calculating. ",date())
  }
  
  bMNTD.rand<-array(unlist(bMNTD.rand),dim=c(nrow(bMNTD.rand[[1]]),ncol(bMNTD.rand[[1]]),length(bMNTD.rand)))
  
  bNTI=(bMNTD.obs-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),stats::sd))
  diag(bNTI)<-0
  bNTI[is.na(bNTI)]=0
  bNTI[is.na(bMNTD.obs)]=NA
  colnames(bNTI)<-rownames(bNTI)<-samp.name
  
  if(RC)
  {
    bMNTD.obsar=array(bMNTD.obs,dim=dim(bMNTD.rand))
    alpha1=apply(bMNTD.obsar==bMNTD.rand,c(1,2),sum)
    alpha=apply(bMNTD.obsar>bMNTD.rand,c(1,2),sum)
    alpha=(alpha+0.5*alpha1)/rand
    rc.res=2*alpha-1
    colnames(rc.res)<-rownames(rc.res)<-samp.name
    result=rc.res
  }else{
    rc.res=NA
    result=bNTI
  }
  gc()
  
  colnames(result)<-rownames(result)<-samp.name
  
  if(output.dtail)
  {
    colnames(bMNTD.obs)<-rownames(bMNTD.obs)<-samp.name
    
    rownames(bMNTD.rand)=rownames(bMNTD.obs)
    colnames(bMNTD.rand)=colnames(bMNTD.obs)
    bMNTD.randm=sapply(1:(dim(bMNTD.rand)[3]),
                       function(i)
                       {
                         (iCAMP::dist.3col(bMNTD.rand[,,i]))[,3]
                       })
    colnames(bMNTD.randm)=paste0("rand",1:ncol(bMNTD.randm))
    bMNTD.randout=data.frame((iCAMP::dist.3col(bMNTD.rand[,,1]))[,1:2,drop=FALSE],
                             bMNTD.randm,stringsAsFactors = FALSE)
    output=list(bNTI=bNTI,RC.bMNTD=rc.res,bMNTD=bMNTD.obs,bMNTD.rand=bMNTD.randout)
  }else{
    output=result
  }
  output
}