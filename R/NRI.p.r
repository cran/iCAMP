NRI.p<-function(comm, dis, nworker=4, memo.size.GB=50,
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
    comm=spc$comm;dis=spc$dis
  }
  ## randomization function ##
  sp.num=ncol(comm)
  
  requireNamespace("permute")
  permat=permute::shuffleSet(sp.num,rand)
  if(nrow(permat)<rand)
  {
    permat=rbind(permat,1:sp.num)
    message("total possible permutation number is ",nrow(permat),", less than random time setting.")
  }
  permat=lapply(1:nrow(permat), function(i,per){per[i,]},per=permat)
  gc()
  
  MPD.random<-function(permati,diss,com,weighted)
  {
    requireNamespace("iCAMP")
    diss.rand=diss
    rand.name=colnames(diss)[permati]
    colnames(diss.rand)=rand.name
    rownames(diss.rand)=rand.name
    MPD.rand<-as.matrix(iCAMP::mpdn(com, diss.rand, abundance.weighted = weighted))
    MPD.rand
  }
  # calculate across all samples #
  if(!silent){message("Now calculating observed MPD. Begin at ", date(),". Please wait...")}
  MPD.obs<-as.matrix(iCAMP::mpdn(comm, dis, abundance.weighted = weighted)) # calculate observed MPD.
  spname=colnames(comm)
  gc()
  if(nworker==1)
  {
    MPD.rand=lapply(permat,MPD.random,diss=dis,com=comm,weighted=weighted)
  }else{
    requireNamespace("parallel")
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
    if(!silent){message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")}
    MPD.rand<-parallel::parLapply(c1,permat,MPD.random,diss=dis,com=comm,weighted=weighted)
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