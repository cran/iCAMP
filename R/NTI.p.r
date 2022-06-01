NTI.p<-function(comm, dis, nworker=4, memo.size.GB=50, 
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
    dis=spc$dis
  }
 
  ## randomization function ##
  
  requireNamespace("permute")
  sp.num=ncol(comm)
  permat=permute::shuffleSet(sp.num,rand)
  if(nrow(permat)<rand)
  {
    permat=rbind(permat,1:sp.num)
    if(!silent){message("total possible permutation number is ",nrow(permat),", less than random time setting.")}
    rand=nrow(permat)
  }
  permat=lapply(1:nrow(permat), function(i,per){per[i,]},per=permat)
  gc()
  
  MNTD.random<-function(permati,diss,com,weighted)
  {
    requireNamespace("iCAMP")
    diss.rand=diss
    rand.name=colnames(diss)[permati]
    colnames(diss.rand)=rand.name
    rownames(diss.rand)=rand.name
    gc()
    MNTD.rand<-as.matrix(iCAMP::mntdn(com, diss.rand, abundance.weighted = weighted))
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
  MNTD.rand<-parallel::parLapply(c1,permat,MNTD.random,diss=dis,com=comm,weighted=weighted)
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