NRI.cm<-function(comm, dis, meta.group=NULL, meta.spool=NULL,
                nworker=4, memo.size.GB=50,
                weighted=c(TRUE,FALSE), check.name=TRUE,
                rand=1000,output.MPD=c(FALSE,TRUE),silent=FALSE,
                sig.index=c("SES","NRI","Confidence","RC","all"))
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
  output.MPD=output.MPD[1]
  if(!(sig.index[1] %in% c("SES","NRI","Confidence","RC","all"))){stop("wrong sig.index for NRI.p.")}
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
  sp.num=ncol(comm)
  
  MPD.random<-function(i,dis,comms,weighted,
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
    MPD.rand<-as.matrix(iCAMP::mpdn(comr, dis, abundance.weighted = weighted))
    MPD.rand
  }
  
  # calculate across all samples #
  if(!silent){message("Now calculating observed MPD. Begin at ", date(),". Please wait...")}
  MPD.obs<-as.matrix(iCAMP::mpdn(comm, dis, abundance.weighted = weighted)) # calculate observed MPD.
  spname=colnames(comm)
  gc()
  if(nworker==1)
  {
    MPD.rand=lapply(1:rand,MPD.random,diss=dis,com=comm,weighted=weighted)
  }else{
    requireNamespace("parallel")
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
    if(!silent){message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")}
    MPD.rand<-parallel::parLapply(c1,1:rand,MPD.random,dis,comms,weighted,
                                  meta.group,meta.spool,meta.spall)
    parallel::stopCluster(c1)
    gc()
  }
  
  MPD.rand<-array(unlist(MPD.rand),dim=c(nrow(MPD.rand[[1]]),ncol(MPD.rand[[1]]),length(MPD.rand)))
  if(output.MPD)
  {
    MPDrandm=matrix(MPD.rand,nrow = dim(MPD.rand)[1],ncol = dim(MPD.rand)[3])
    rownames(MPDrandm)=rownames(MPD.obs)
    colnames(MPDrandm)=paste0("rand",1:ncol(MPDrandm))
  }
  
  if(sig.index[1] %in% c("SES","NRI","all"))
  {
    NRI=(apply(MPD.rand,c(1,2),mean)-MPD.obs)/(apply(MPD.rand,c(1,2),stats::sd))
    NRI[is.na(NRI)]=0
  }
  
  if(sig.index[1] %in% c("Confidence","all"))
  {
    MPD.obsar=array(MPD.obs,dim=dim(MPD.rand))
    alpha=(apply(MPD.obsar>MPD.rand,c(1,2),sum))/rand
    alpha2=(apply(MPD.obsar<MPD.rand,c(1,2),sum))/rand
    alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
    alpha[is.na(alpha)]=0
    rownames(alpha)=rownames(MPD.obs)
  }
  
  if(sig.index[1] %in% c("RC","all"))
  {
    MPD.obsar=array(MPD.obs,dim=dim(MPD.rand))
    alphax=(apply(MPD.obsar>MPD.rand,c(1,2),sum))/rand
    alphax2=(apply(MPD.obsar==MPD.rand,c(1,2),sum))/rand
    alphax=(alphax+0.5*alphax2)
    RC=2*alphax-1
  }
  
  if(sig.index[1] %in% c("SES","NRI"))
  {
    if(output.MPD)
    {
      output=list(index=data.frame(NRI=NRI,stringsAsFactors = FALSE),
                  MPD.obs=data.frame(MPD=MPD.obs),MPD.rand=MPDrandm)
    }else{
      output=data.frame(NRI=NRI,stringsAsFactors = FALSE)
    }
  }else if(sig.index[1]=="Confidence"){
    if(output.MPD)
    {
      output=list(index=data.frame(CMPD=alpha,stringsAsFactors = FALSE),
                  MPD.obs=data.frame(MPD=MPD.obs),MPD.rand=MPDrandm)
    }else{output=data.frame(CMPD=alpha,stringsAsFactors = FALSE)}
  }else if(sig.index[1]=="RC"){
    if(output.MPD)
    {
      output=list(index=data.frame(RCMPD=RC,stringsAsFactors = FALSE),
                  MPD.obs=data.frame(MPD=MPD.obs),MPD.rand=MPDrandm)
    }else{output=data.frame(RCMPD=RC,stringsAsFactors = FALSE)}
  }else if(sig.index[1]=="all"){
    if(output.MPD)
    {
      output=list(SES=data.frame(NRI=NRI,stringsAsFactors = FALSE),
                  Confidence=data.frame(CMPD=alpha,stringsAsFactors = FALSE),
                  RC=data.frame(RCMPD=RC,stringsAsFactors = FALSE),
                  MPD.obs=data.frame(MPD=MPD.obs),MPD.rand=MPDrandm)
    }else{
      output=list(SES=data.frame(NRI=NRI,stringsAsFactors = FALSE),
                  Confidence=data.frame(CMPD=alpha,stringsAsFactors = FALSE),
                  RC=data.frame(RCMPD=RC,stringsAsFactors = FALSE))
    }
  }
  output
}