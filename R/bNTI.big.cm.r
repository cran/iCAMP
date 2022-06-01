bNTI.big.cm<-function(comm, meta.group=NULL, meta.spool=NULL, pd.desc="pd.desc",
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
  
  # check sample and species IDs. Split comm into comms of different metacommunities.
  if(!is.null(meta.group))
  {
    sampc=iCAMP::match.name(rn.list = list(comm=comm,meta.group=meta.group))
    comm=sampc$comm
    meta.group=sampc$meta.group
    meta.lev=unique(meta.group[,1])
  }else{
    meta.group=data.frame(metagrp=rep("Meta",nrow(comm)),stringsAsFactors = FALSE)
    rownames(meta.group)=rownames(comm)
    meta.lev="Meta"
  }
  if(spname.check)
  {
    #aslist<-function(a){if(is.null(a)){NULL}else{out=list(a);names(out)=deparse(substitute(a));out}}
    #sampc=iCAMP::match.name(rn.list = c(aslist(comm),aslist(meta.group)))
    #comm=sampc$comm
    #if(!is.null(meta.group)) meta.group=sampc$meta.group
    check.sp=match.name(name.check = pd.spname,cn.list = list(comm=comm))
    comm=check.sp$comm
  }
  if(!is.null(meta.spool))
  {
    if(!is.list(meta.spool))
    {
      if(length(meta.lev)>1){stop('meta.spool needs to specify species pool for each metacommunity in meta.group.')}
      meta.spool=list(Meta=meta.spool)
    }
    if(sum(!(unique(unlist(meta.spool)) %in% pd.spname))>0){stop('meta.spool has some species not included in pd.spname.')}
    if(sum(!(meta.lev %in% names(meta.spool)))>0)
    {
      stop('meta.spool element names must be the same as metacommunity names in meta.group.')
    }else{
      meta.spool=meta.spool[match(meta.lev,names(meta.spool))]
    }
    comms=lapply(meta.lev,function(mi){sampi=rownames(meta.group)[which(meta.group[,1]==mi)];comi=comm[which(rownames(comm) %in% sampi),,drop=FALSE];comi[,colSums(comi)>0,drop=FALSE]})
    names(comms)=meta.lev
    msplen=sapply(meta.spool,length)
    mgplen=sapply(comms,ncol)
    if(sum((msplen-mgplen)<0)>0){stop('meta.spool shouldnot have less species than observed species in a metacommunity.')}
  }else{
    meta.spool=list(Meta=colnames(comm))
    comms=list(Meta=comm)
  }
  meta.spall=unique(unlist(meta.spool))
  
  samp.name=rownames(comm)
  if(trace){trace.seq=floor(seq(from=1,to=rand,by=50))}else{trace.seq=NULL}
  
  ## randomization function ##
  sp.num=length(pd.spname)
  
  bMNTD.random<-function(i,pd.desc,comms,meta.group,meta.spool,meta.spall,pd.spname,
                         weighted,exclude.consp,pd.wd,trace.seq)
  {
    if(!is.null(trace.seq)){if(i %in% trace.seq) message("Now randomizing i=",i,". ",date())}
    
    comr=matrix(0,nrow=nrow(meta.group),ncol=length(meta.spall))
    rownames(comr)=rownames(meta.group)
    colnames(comr)=meta.spall
    for(j in 1:length(meta.spool))
    {
      comj=comms[[j]]
      if(ncol(comj)>0)
      {
        spj.rd=sample(meta.spool[[j]],ncol(comj))
        comr[match(rownames(comj),rownames(comr)),match(spj.rd,colnames(comr))]=comj
      }
    }
    comr=comr[,colSums(comr)>0,drop=FALSE]
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
    bMNTD.rand=lapply(1:rand,bMNTD.random,
                      pd.desc=pd.desc,comms=comms,meta.group=meta.group,
                      meta.spool=meta.spool,meta.spall=meta.spall,pd.spname=pd.spname,
                      weighted=weighted,exclude.consp=exclude.consp,
                      pd.wd=pd.wd,trace.seq=trace.seq)
  }else{
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
    if(trace) message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
    bMNTD.rand<-parallel::parLapply(c1,1:rand,bMNTD.random,
                                    pd.desc=pd.desc,comms=comms,meta.group=meta.group,
                                    meta.spool=meta.spool,meta.spall=meta.spall,pd.spname=pd.spname,
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