bNTIn.p<-function(comm, dis, nworker=4, memo.size.GB=50,
                  weighted=c(TRUE,FALSE),exclude.consp=FALSE,rand=1000,
                  output.bMNTD=c(FALSE,TRUE),sig.index=c("SES","Confidence","RC","bNTI"),
                  unit.sum=NULL,correct.special=FALSE,detail.null=FALSE,
                  special.method=c("MNTD", "MPD", "both"),
                  ses.cut=1.96,rc.cut=0.95,conf.cut=0.975)
{
  #v20200725 add conf.cut, detail.null. change RC to sig.index.
  #load package
  requireNamespace("parallel")
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memo.size.GB*1024)
    {
      memotry=try(utils::memory.limit(size=memo.size.GB*1024),silent = TRUE)
      if(class(memotry)=="try-error"){warning(memotry[1])}
    }
  }
  samp.name=rownames(comm)
  weighted=weighted[1]
  output.bMNTD=output.bMNTD[1]
  sig.index=sig.index[1]
  special.method=special.method[1]
  if(!(sig.index %in% c("SES","Confidence","RC","bNTI"))){stop("wrong sig.index for bNTIn.p.")}
  
  ## randomization function ##
  bMNTD.random<-function(i,diss,com,weighted,exclude.consp,unit.sum)
  {
    requireNamespace("iCAMP")
    diss.rand=diss
    rand.name=sample(colnames(diss))
    colnames(diss.rand)=rand.name
    rownames(diss.rand)=rand.name
    gc()
    bMNTD.rand<-as.matrix(iCAMP::bmntd(com, diss.rand, abundance.weighted = weighted, exclude.conspecifics = exclude.consp,unit.sum=unit.sum))
    bMNTD.rand
  }
  
  
  # calculate across all samples #
  message("Now calculating observed betaMNTD. Begin at ", date(),". Please wait...")
  gc()
  bMNTD.obs<-as.matrix(iCAMP::bmntd(comm, dis, abundance.weighted = weighted, exclude.conspecifics = exclude.consp,unit.sum=unit.sum)) # calculate observed betaMNTD.
  spname=colnames(comm)
  gc()
  c1<-parallel::makeCluster(nworker,type="PSOCK")
  message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
  bMNTD.rand<-parallel::parLapply(c1,1:rand,bMNTD.random,diss=dis,com=comm,weighted=weighted,exclude.consp=exclude.consp,unit.sum=unit.sum)
  parallel::stopCluster(c1)
  gc()
  bMNTD.rand<-array(unlist(bMNTD.rand),dim=c(nrow(bMNTD.rand[[1]]),ncol(bMNTD.rand[[1]]),length(bMNTD.rand)))
  
  if(sig.index=="RC")
  {
    bMNTD.obsar=array(bMNTD.obs,dim=dim(bMNTD.rand))
    alpha1=apply(bMNTD.obsar==bMNTD.rand,c(1,2),sum)
    alpha=apply(bMNTD.obsar>bMNTD.rand,c(1,2),sum)
    alpha=(alpha+0.5*alpha1)/rand
    result=2*alpha-1
  }else if(sig.index %in% c("SES","bNTI")){
    bNTI=(bMNTD.obs-apply(bMNTD.rand,c(1,2),mean))/(apply(bMNTD.rand,c(1,2),stats::sd))
    diag(bNTI)<-0
    result=bNTI
  }else if(sig.index=="Confidence"){
    bMNTD.obsar=array(bMNTD.obs,dim=dim(bMNTD.rand))
    alpha=(apply(bMNTD.obsar>bMNTD.rand,c(1,2),sum))/rand
    alpha2=(apply(bMNTD.obsar<bMNTD.rand,c(1,2),sum))/rand
    alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
    result=alpha
  }
  colnames(result)<-rownames(result)<-samp.name
  gc()
  
  if(correct.special)
  {
    # some special cases
    message("Now fixing special cases. Begin at ",date(),". Please wait...")
    sdm=(apply(bMNTD.rand,c(1,2),stats::sd))
    diag(sdm)<-NA
    if(detail.null){special.ses=result;special.ses[]=0;special.rc<-special.conf<-special.ses}
    if(sum(sdm==0,na.rm = TRUE)>0)
    {
      rownames(sdm)<-colnames(sdm)<-rownames(bMNTD.obs)
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
                       rand = rand,output.MPD = FALSE,sig.index=sig.nti)
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
      
      if(length(samp1.id)>0 | length(samp.sg)>0)
      {
        rcm=(iCAMP::RC.pc(comm=comm,rand=rand,na.zero=TRUE,nworker=nworker,
                          memory.G=memo.size.GB,weighted=weighted,
                          unit.sum=unit.sum,silent=TRUE))$index
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
  if(output.bMNTD[1])
  {
    colnames(bMNTD.obs)<-rownames(bMNTD.obs)<-samp.name
    output=c(output,list(betaMNTD.obs=bMNTD.obs))
  }
  if(detail.null)
  {
    rownames(bMNTD.rand)=rownames(bMNTD.obs)
    colnames(bMNTD.rand)=colnames(bMNTD.obs)
    bMNTD.randm=sapply(1:(dim(bMNTD.rand)[3]),
                       function(i)
                       {
                         (iCAMP::dist.3col(bMNTD.rand[,,i]))[,3]
                       })
    colnames(bMNTD.randm)=paste0("rand",1:ncol(bMNTD.randm))
    if(correct.special)
    {
      special.crct=list(special.ses=special.ses,special.rc=special.rc,special.conf=special.conf)
    }else{special.crct=NULL}
    output=c(output,list(rand=data.frame((iCAMP::dist.3col(bMNTD.rand[,,1]))[,1:2,drop=FALSE],
                                         bMNTD.randm,stringsAsFactors = FALSE),
                         special.crct=special.crct))
  }
  output
}
