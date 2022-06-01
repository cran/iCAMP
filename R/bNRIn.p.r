bNRIn.p<-function(comm, dis, nworker=4, memo.size.GB=50,
                  weighted=c(TRUE,FALSE),rand=1000,output.bMPD=c(FALSE,TRUE),
                  sig.index=c("SES","Confidence","RC","bNRI"),unit.sum=NULL,
                  correct.special=FALSE,detail.null=FALSE,
                  special.method=c("MPD","MNTD","both"),
                  ses.cut=1.96,rc.cut=0.95,conf.cut=0.975,
                  dirichlet = FALSE)
{
  #v20200727 add conf.cut, detail.null. change RC to sig.index.
  #load package
  requireNamespace("parallel")
  requireNamespace("permute")
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memo.size.GB*1024)
    {
      memotry=try(utils::memory.limit(size=memo.size.GB*1024),silent = TRUE)
      if(inherits(memotry,"try-error")){warning(memotry[1])}
    }
  }
  
  
  weighted=weighted[1]
  output.bMPD=output.bMPD[1]
  sig.index=sig.index[1]
  special.method=special.method[1]
  if(!(sig.index %in% c("SES","Confidence","RC","bNRI"))){stop("wrong sig.index for bNRIn.p.")}
  
  # match
  if(sum(colnames(comm)!=rownames(dis))>0)
  {
    sp.name=intersect(colnames(comm),rownames(dis))
    comm=comm[,match(sp.name,colnames(comm))]
    dis=dis[match(sp.name,rownames(dis)),match(sp.name,rownames(dis))]
  }
  #print(dim(comm))
  #print(dim(dis))
  samp.name=rownames(comm)
  
  sp.name=rownames(dis)
  sp.num=nrow(dis)
  ## randomization function ##
  perm=permute::shuffleSet(sp.num,rand)
  if(sp.num==2){perm=matrix(c(2,1),nrow=1)}
  if(nrow(perm)<rand)
  {
    perm=rbind(perm,1:sp.num)
    rand=nrow(perm)
  }
  
  bMPD.random<-function(i,diss,com,weighted,perm,unit.sum)
  {
    requireNamespace("iCAMP")
    diss.rand=diss[perm[i,],perm[i,]]
    colnames(diss.rand)<-rownames(diss.rand)<-rownames(diss)
    gc()
    bMPD.rand<-as.matrix(iCAMP::bmpd(com, diss.rand, abundance.weighted = weighted, unit.sum=unit.sum))
    bMPD.rand
  }
  
  
  # calculate across all samples #
  message("Now calculating observed betaMPD. Begin at ", date(),". Please wait...")
  gc()
  bMPD.obs<-as.matrix(iCAMP::bmpd(comm, dis, abundance.weighted = weighted, unit.sum=unit.sum)) # calculate observed betaMPD.
  spname=colnames(comm)
  gc()
  
  c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
  if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
  if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
  message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
  bMPD.rand<-parallel::parLapply(c1,1:rand,bMPD.random,diss=dis,com=comm,weighted=weighted,perm=perm,unit.sum=unit.sum)
  parallel::stopCluster(c1)
  gc()
  bMPD.rand<-array(unlist(bMPD.rand),dim=c(nrow(bMPD.rand[[1]]),ncol(bMPD.rand[[1]]),length(bMPD.rand)))
  if(sig.index=="RC")
  {
    bMPD.obsar=array(bMPD.obs,dim=dim(bMPD.rand))
    alpha1=apply(bMPD.obsar==bMPD.rand,c(1,2),sum)
    alpha=apply(bMPD.obsar>bMPD.rand,c(1,2),sum)
    alpha=(alpha+0.5*alpha1)/rand
    result=2*alpha-1
  }else if(sig.index %in% c("SES","bNRI")){
    bNRI=(bMPD.obs-apply(bMPD.rand,c(1,2),mean))/(apply(bMPD.rand,c(1,2),stats::sd))
    result=bNRI
    gc()
    diag(result)<-0
  }else if(sig.index=="Confidence"){
    bMPD.obsar=array(bMPD.obs,dim=dim(bMPD.rand))
    alpha=(apply(bMPD.obsar>bMPD.rand,c(1,2),sum))/rand
    alpha2=(apply(bMPD.obsar<bMPD.rand,c(1,2),sum))/rand
    alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
    result=alpha
  }
  colnames(result)<-rownames(result)<-samp.name
  
  if(correct.special)
  {
    # some special cases
    message("Now fixing special cases. Begin at ",date(),". Please wait...")
    sdm=(apply(bMPD.rand,c(1,2),stats::sd))
    diag(sdm)<-NA
    if(detail.null){special.ses=result;special.ses[]=0;special.rc<-special.conf<-special.ses}
    if(sum(sdm==0,na.rm = TRUE)>0)
    {
      rownames(sdm)<-colnames(sdm)<-rownames(bMPD.obs)
      sdc=dist.3col(sdm)
      sdc0=sdc[which(sdc[,3]==0),,drop=FALSE]
      samp.sd0=unique(as.vector(as.matrix(sdc0[,1:2])))
      samp0=rownames(comm)[which(rowSums(comm)==0)]
      samp.ck=setdiff(samp.sd0,samp0)
      if(length(samp.ck)>0)
      {
        comm.ck=comm[which(rownames(comm) %in% samp.ck),,drop=FALSE]
        samp.sg=character(0)
        if(detail.null)
        {
          samp.sg.ses<-samp.sg.rc<-samp.sg.conf<-character(0)
          sig.nti="all"
        }else{
          if(sig.index %in% c("bNRI","SES")){sig.nti="SES";nti.cut=ses.cut}else if(sig.index=="RC"){sig.nti="RC";nti.cut=rc.cut}else if(sig.index=="Confidence"){sig.nti="Confidence";nti.cut=conf.cut}
        }
        if(special.method[1] %in% c("MNTD","NTI","both"))
        {
          nti.ck=iCAMP::NTI.p(comm = comm.ck, dis = dis, nworker = nworker,
                       memo.size.GB = memo.size.GB,weighted = weighted,
                       rand = rand,output.MNTD = FALSE,sig.index=sig.nti)
          if(detail.null)
          {
            samp.sg.ses=c(samp.sg.ses,rownames(nti.ck$SES)[which(abs(nti.ck$SES[,1])>ses.cut)])
            samp.sg.rc=c(samp.sg.rc,rownames(nti.ck$RC)[which(abs(nti.ck$RC[,1])>rc.cut)])
            samp.sg.conf=c(samp.sg.rc,rownames(nti.ck$Confidence)[which(abs(nti.ck$Confidence[,1])>conf.cut)])
          }else{
            samp.sg=c(samp.sg,rownames(nti.ck)[which(abs(nti.ck[,1])>nti.cut)])
          }
        }
        
        if(special.method[1] %in% c("MPD","NRI","both"))
        {
          nri.ck=iCAMP::NRI.p(comm = comm.ck, dis = dis, nworker = nworker,
                       memo.size.GB = memo.size.GB,weighted = weighted,
                       rand = rand,output.MPD = FALSE,sig.index = sig.nti)
          if(detail.null)
          {
            samp.sg.ses=c(samp.sg.ses,rownames(nri.ck$SES)[which(abs(nri.ck$SES[,1])>ses.cut)])
            samp.sg.rc=c(samp.sg.rc,rownames(nri.ck$RC)[which(abs(nri.ck$RC[,1])>rc.cut)])
            samp.sg.conf=c(samp.sg.rc,rownames(nri.ck$Confidence)[which(abs(nri.ck$Confidence[,1])>conf.cut)])
          }else{
            samp.sg=c(samp.sg,rownames(nri.ck)[which(abs(nri.ck[,1])>nti.cut)])
          }
        }
        
        if(detail.null)
        {
          if(sig.index %in% c("bNRI","SES")){samp.sg=samp.sg.ses}else if(sig.index=="RC"){samp.sg=samp.sg.rc}else if(sig.index=="Confidence"){samp.sg=samp.sg.conf}
          samp.sg.sum=length(samp.sg.ses)+length(samp.sg.rc)+length(samp.sg.conf)
        }else{samp.sg.sum=0}
        
      }else{samp.sg=character(0);samp.sg.sum=0}
      
      samp1.id=which(rowSums(comm>0)==1)
      
      if(length(samp1.id)>0 | length(samp.sg)>0 | samp.sg.sum>0)
      {
        rcm=(iCAMP::RC.pc(comm=comm,rand=rand,na.zero=TRUE,nworker=nworker,
                          memory.G=memo.size.GB,weighted=weighted,
                          unit.sum=unit.sum,silent=TRUE,dirichlet = dirichlet))$index
        
        if(length(samp1.id)>0)
        {
          requireNamespace("vegan")
          BCm=as.matrix(vegan::vegdist(comm,method = "euclidean",binary = TRUE))
          diag(BCm)=NA
          id.samesp=which(BCm==0,arr.ind = TRUE)
          id.same1=id.samesp[which((id.samesp[,1] %in% samp1.id)&(id.samesp[,2] %in% samp1.id)),,drop=FALSE]
          
          if(sig.index %in% c("RC","Confidence"))
          {
            result[id.same1]=(1-(2*(rcm[id.same1]<=0)))*1.1
          }else if(sig.index %in% c("SES","bNTI")){
            result[id.same1]=(1-(2*(rcm[id.same1]<=0)))*99
          }
          
          if(detail.null)
          {
            special.ses[id.same1]=(1-(2*(rcm[id.same1]<=0)))*99
            special.rc[id.same1]<-special.conf[id.same1]<-((1-(2*(rcm[id.same1]<=0)))*1.1)
          }
        }
        sd0.rcn=which(sdm==0,arr.ind = TRUE)
        if(length(samp.sg)>0)
        {
          sd0.sg=sd0.rcn[which((rownames(sdm)[sd0.rcn[,1]] %in% samp.sg)|(rownames(sdm)[sd0.rcn[,2]] %in% samp.sg)),,drop=FALSE]
          if(sig.index %in% c("RC","Confidence"))
          {
            result[sd0.sg]=(1-(2*(rcm[sd0.sg]<=0)))*1.1
          }else if(sig.index %in% c("SES","bNTI")){
            result[sd0.sg]=(1-(2*(rcm[sd0.sg]<=0)))*99
          }
        }
        
        if(detail.null & samp.sg.sum>0)
        {
          sd0.sg.ses=sd0.rcn[which((rownames(sdm)[sd0.rcn[,1]] %in% samp.sg.ses)|(rownames(sdm)[sd0.rcn[,2]] %in% samp.sg.ses)),,drop=FALSE]
          special.ses[sd0.sg.ses]=(1-(2*(rcm[sd0.sg.ses]<=0)))*99
          sd0.sg.rc=sd0.rcn[which((rownames(sdm)[sd0.rcn[,1]] %in% samp.sg.rc)|(rownames(sdm)[sd0.rcn[,2]] %in% samp.sg.rc)),,drop=FALSE]
          special.rc[sd0.sg.rc]=(1-(2*(rcm[sd0.sg.rc]<=0)))*1.1
          sd0.sg.conf=sd0.rcn[which((rownames(sdm)[sd0.rcn[,1]] %in% samp.sg.conf)|(rownames(sdm)[sd0.rcn[,2]] %in% samp.sg.conf)),,drop=FALSE]
          special.conf[sd0.sg.conf]=(1-(2*(rcm[sd0.sg.conf]<=0)))*1.1
        }
      }
    }
  }
  
  result[is.na(result)]=0
  
  output=list(index=result,sig.index=sig.index)
  if(output.bMPD[1])
  {
    colnames(bMPD.obs)<-rownames(bMPD.obs)<-samp.name
    output=c(output,list(betaMPD.obs=bMPD.obs))
  }
  if(detail.null)
  {
    rownames(bMPD.rand)=rownames(bMPD.obs)
    colnames(bMPD.rand)=colnames(bMPD.obs)
    bMPD.randm=matrix(sapply(1:(dim(bMPD.rand)[3]),
                       function(i)
                       {
                         (iCAMP::dist.3col(bMPD.rand[,,i]))[,3]
                       }),ncol=(dim(bMPD.rand)[3]))
    colnames(bMPD.randm)=paste0("rand",1:ncol(bMPD.randm))
    if(correct.special)
    {
      special.crct=list(special.ses=special.ses,special.rc=special.rc,special.conf=special.conf)
    }else{special.crct=NULL}
    output=c(output,list(rand=data.frame((iCAMP::dist.3col(bMPD.rand[,,1]))[,1:2,drop=FALSE],
                                         bMPD.randm,stringsAsFactors = FALSE),
                         special.crct=special.crct))
  }
  output
}
