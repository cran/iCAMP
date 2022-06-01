NTI.cm<-function(comm, dis, meta.group=NULL,meta.spool=NULL,
                nworker=4, memo.size.GB=50,
                weighted=c(TRUE,FALSE),rand=1000,
                check.name=TRUE,output.MNTD=c(FALSE,TRUE),
                sig.index=c("SES","NTI","Confidence","RC","all"),silent=FALSE)
{
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memo.size.GB*1024)
    {
      memotry=try(utils::memory.limit(size=memo.size.GB*1024),silent = TRUE)
      if(inherits(memotry,"try-error")){warning(memotry[1])}
    }
  }
  weighted=weighted[1]
  output.MNTD=output.MNTD[1]
  if(!(sig.index[1] %in% c("SES","NTI","Confidence","RC","all"))){stop("wrong sig.index for NTI.p.")}
  # match
  if(check.name)
  {
    spc=iCAMP::match.name(cn.list=list(comm=comm),both.list=list(dis=dis))
    comm=spc$comm
    #dis=spc$dis
  }
  
  ####################
  ## match IDs
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
  comms=lapply(meta.lev,function(mi){sampi=rownames(meta.group)[which(meta.group[,1]==mi)];comi=comm[which(rownames(comm) %in% sampi),,drop=FALSE];comi[,colSums(comi)>0,drop=FALSE]})
  names(comms)=meta.lev
  
  if(!is.null(meta.spool))
  {
    if(!is.list(meta.spool))
    {
      if(length(meta.lev)>1){stop('meta.spool needs to specify species pool for each metacommunity in meta.group.')}
      meta.spool=list(Meta=meta.spool)
    }
    if(sum(!(unique(unlist(meta.spool)) %in% rownames(dis)))>0){stop('meta.spool has some species not included in dis.')}
    if(sum(!(meta.lev %in% names(meta.spool)))>0)
    {
      stop('meta.spool element names must be the same as metacommunity names in meta.group.')
    }else{
      meta.spool=meta.spool[match(meta.lev,names(meta.spool))]
    }
    msplen=sapply(meta.spool,length)
    mgplen=sapply(comms,ncol)
    if(sum((msplen-mgplen)<0)>0){stop('meta.spool shouldnot have less species than observed species in a metacommunity.')}
  }else{
    meta.spool=lapply(1:length(comms),function(i){colnames(comms[[i]])})
    names(meta.spool)=names(comms)
  }
  
  meta.spall=unique(unlist(meta.spool))
  
  ## randomization function ##
  
  requireNamespace("permute")
  sp.num=ncol(comm)
  gc()
  
  MNTD.random<-function(i,dis,comms,weighted,
                        meta.group,meta.spool,meta.spall)
  {
    requireNamespace("iCAMP")
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
    MNTD.rand<-as.matrix(iCAMP::mntdn(comr, dis, abundance.weighted = weighted))
    MNTD.rand
  }

  # calculate across all samples #
  if(!silent){message("Now calculating observed MNTD. Begin at ", date(),". Please wait...")}
  MNTD.obs<-as.matrix(iCAMP::mntdn(comm=comm, pd=dis, abundance.weighted = weighted)) # calculate observed MNTD.
  gc()
  
  requireNamespace("parallel")
  c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
  if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
  if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
  if(!silent){message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")}
  MNTD.rand<-parallel::parLapply(c1,1:rand,MNTD.random,dis,comms,weighted,
                                 meta.group,meta.spool,meta.spall)
  parallel::stopCluster(c1)
  MNTD.rand<-array(unlist(MNTD.rand),dim=c(nrow(MNTD.rand[[1]]),ncol(MNTD.rand[[1]]),length(MNTD.rand)))
  gc()
  rownames(MNTD.rand)=rownames(MNTD.obs)
  
  if(output.MNTD)
  {
    MNTDrandm=matrix(MNTD.rand,nrow = dim(MNTD.rand)[1],ncol = dim(MNTD.rand)[3])
    rownames(MNTDrandm)=rownames(MNTD.obs)
    colnames(MNTDrandm)=paste0("rand",1:ncol(MNTDrandm))
  }
  
  if(sig.index[1] %in% c("SES","NTI", "all"))
  {
    NTI=(apply(MNTD.rand,c(1,2),mean)-MNTD.obs)/(apply(MNTD.rand,c(1,2),stats::sd))
    NTI[is.na(NTI)]=0
  }
  
  if(sig.index[1] %in% c("Confidence","all"))
  {
    MNTD.obsar=array(MNTD.obs,dim=dim(MNTD.rand))
    alpha=(apply(MNTD.obsar>MNTD.rand,c(1,2),sum))/rand
    alpha2=(apply(MNTD.obsar<MNTD.rand,c(1,2),sum))/rand
    alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
    alpha[is.na(alpha)]=0
    rownames(alpha)=rownames(MNTD.obs)
  }
  
  if(sig.index[1] %in% c("RC","all"))
  {
    MNTD.obsar=array(MNTD.obs,dim=dim(MNTD.rand))
    alphax=(apply(MNTD.obsar>MNTD.rand,c(1,2),sum))/rand
    alphax2=(apply(MNTD.obsar==MNTD.rand,c(1,2),sum))/rand
    alphax=(alphax+0.5*alphax2)
    RC=2*alphax-1
  }
  
  if(sig.index[1] %in% c("SES","NTI"))
  {
    if(output.MNTD)
    {
      output=list(index=data.frame(NTI=NTI,stringsAsFactors = FALSE),
                  MNTD.obs=data.frame(MNTD=MNTD.obs),MNTD.rand=MNTDrandm)
    }else{output=data.frame(NTI=NTI,stringsAsFactors = FALSE)}
  }else if(sig.index[1] == "Confidence"){
    if(output.MNTD)
    {
      output=list(index=data.frame(CMNTD=alpha,stringsAsFactors = FALSE),
                  MNTD.obs=data.frame(MNTD=MNTD.obs),MNTD.rand=MNTDrandm)
    }else{output=data.frame(CMNTD=alpha,stringsAsFactors = FALSE)}
  }else if(sig.index[1] == "RC"){
    if(output.MNTD)
    {
      output=list(index=data.frame(RCMNTD=RC,stringsAsFactors = FALSE),
                  MNTD.obs=data.frame(MNTD=MNTD.obs),MNTD.rand=MNTDrandm)
    }else{output=data.frame(RCMNTD=RC,stringsAsFactors = FALSE)}
  }else if(sig.index[1] == "all"){
    if(output.MNTD)
    {
      output=list(SES=data.frame(NTI=NTI,stringsAsFactors = FALSE),
                  Confidence=data.frame(CMNTD=alpha,stringsAsFactors = FALSE),
                  RC=data.frame(RCMNTD=RC,stringsAsFactors = FALSE),
                  MNTD.obs=data.frame(MNTD=MNTD.obs),MNTD.rand=MNTDrandm)
    }else{
      output=list(SES=data.frame(NTI=NTI,stringsAsFactors = FALSE),
                  Confidence=data.frame(CMNTD=alpha,stringsAsFactors = FALSE),
                  RC=data.frame(RCMNTD=RC,stringsAsFactors = FALSE))
    }
  }
  output
}