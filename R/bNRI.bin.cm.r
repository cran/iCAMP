bNRI.bin.cm<-function(comm, meta.group=NULL,meta.spool=NULL,
                      meta.frequency=NULL,meta.ab=NULL,
                      pd.desc, pd.spname, pd.wd, pdid.bin, sp.bin, 
                      spname.check=FALSE,nworker=4, memo.size.GB=50, weighted=c(TRUE,FALSE),
                      rand=1000,output.bMPD=FALSE,sig.index=c("SES","Confidence","RC","bNRI"),
                      unit.sum=NULL,correct.special=FALSE,detail.null=FALSE,
                      special.method=c("MPD","MNTD","both"),ses.cut=1.96,rc.cut=0.95,conf.cut=0.975,
                      dirichlet = FALSE)
{
  #v20200727 add conf.cut, detail.null. change RC to sig.index.
  #load package
  requireNamespace("parallel")
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memo.size.GB*1024)
    {
      memotry=try(utils::memory.limit(size=memo.size.GB*1024),silent = TRUE)
      if(inherits(memotry,"try-error")){warning(memotry[1])}
    }
  }
  
  ###############################
  # match
  if(spname.check)
  {
    sp.check=iCAMP::match.name(name.check=pd.spname,cn.list = list(comm=comm))
    comm=sp.check$comm
  }
  samp.name=rownames(comm)
  sp.name=colnames(comm)
  sp.num=ncol(comm)
  
  ######################
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
    if(sum(!(unique(unlist(meta.spool)) %in% pd.spname))>0){stop('meta.spool has some species not included in pd.spname.')}
    if(sum(!(meta.lev %in% names(meta.spool)))>0)
    {
      stop('meta.spool element names must be the same as metacommunity names in meta.group.')
    }else{
      meta.spool=meta.spool[match(meta.lev,names(meta.spool))]
    }
    msplen=sapply(meta.spool,length)
    mgplen=sapply(comms,ncol)
    if(sum((msplen-mgplen)<0)>0){stop('meta.spool shouldnot have less species than observed species in a metacommunity.')}
  }
  
  if(correct.special)
  {
    if(!is.null(meta.frequency))
    {
      if(sum(!(colnames(comm) %in% colnames(meta.frequency)))){stop('comm has some species not included in meta.frequence.')}
      if(sum(!(colnames(meta.frequency) %in% pd.spname))>0){stop('meta.frequency has some species not included in pd.spname.')}
      if(sum(!(meta.lev %in% rownames(meta.frequency)))>0)
      {
        stop('meta.frequency rownames must be the same as metacommunity names in meta.group.')
      }else{
        meta.frequency=meta.frequency[match(meta.lev,rownames(meta.frequency)),,drop=FALSE]
      }
    }else{
      meta.frequency=matrix(0,nrow = length(meta.lev),ncol = ncol(comm))
      rownames(meta.frequency)=meta.lev
      colnames(meta.frequency)=colnames(comm)
      for(i in 1:length(comms))
      {
        meta.frequency[i,match(colnames(comms[[i]]),colnames(comm))]=colSums(comms[[i]]>0)
      }
    }
    
    if(is.null(meta.spool))
    {
      meta.spool=lapply(1:nrow(meta.frequency),function(i){colnames(meta.frequency)[which(meta.frequency[i,]>0)]})
    }else{
      if(sum(sapply(1:nrow(meta.frequency),function(i){sum(!(colnames(meta.frequency)[which(meta.frequency[i,]>0)] %in% meta.spool[[i]]))}))>0)
      {
        stop('meta.frequency is not consistent with meta.spool')
      }
    }
    # 
    if(!is.null(meta.ab))
    {
      if(sum(!(colnames(comm) %in% colnames(meta.ab)))){stop('comm has some species not included in meta.ab')}
      if(sum(!(colnames(meta.ab) %in% pd.spname))>0){stop('meta.ab has some species not included in pd.spname.')}
      if(sum(!(colnames(meta.frequency) %in% colnames(meta.ab)))){stop('meta.frequency has some species not included in meta.ab')}
      if(sum(!(meta.lev %in% rownames(meta.ab)))>0)
      {
        stop('meta.ab rownames must be the same as metacommunity names in meta.group.')
      }else{
        meta.ab=meta.ab[match(meta.lev,rownames(meta.ab)),match(colnames(meta.frequency),colnames(meta.ab)),drop=FALSE]
      }
    }else{
      meta.ab=matrix(0,nrow = length(meta.lev),ncol = ncol(comm))
      rownames(meta.ab)=meta.lev
      colnames(meta.ab)=colnames(comm)
      for(i in 1:length(comms))
      {
        meta.ab[i,match(colnames(comms[[i]]),colnames(comm))]=colMeans(comms[[i]]/rowSums(comms[[i]]))
      }
    }
  }
  
  if(is.null(meta.spool))
  {
    meta.spool=lapply(1:length(comms),function(i){colnames(comms[[i]])})
    names(meta.spool)=names(comms)
  }
  
  meta.spall=unique(unlist(meta.spool))
  
  ###########################
  
  weighted=weighted[1];output.bMPD=output.bMPD[1]
  sig.index=sig.index[1]
  special.method=special.method[1]
  if(!(sig.index %in% c("SES","Confidence","RC","bNRI"))){stop("wrong sig.index for bNRI.bin.big.")}
  requireNamespace("bigmemory")
  pd=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
  
  
  
  # bin comm and pdist
  bin.lev=levels(as.factor(sp.bin[,1]))
  bin.num=length(bin.lev)
  com.bin=lapply(1:bin.num, function(i){comm[,match(rownames(sp.bin)[which(sp.bin==i)],colnames(comm))]})

  ## randomization function ##
  
  bMPD.random<-function(i,pd.desc, pd.spname, pd.wd, pdid.bin,
                        comms,sp.bin,bin.num,bin.lev,weighted,unit.sum,
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
    comr.bin=lapply(1:bin.num, function(j){comr[,which(colnames(comr) %in% rownames(sp.bin)[which(sp.bin==bin.lev[j])]),drop=FALSE]})
    
    requireNamespace("bigmemory")
    pd=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
    bMPD.rand<-lapply(1:bin.num,
                      function(u)
                      {
                        pd.binu=pd[pdid.bin[[u]],pdid.bin[[u]]]
                        gc()
                        rownames(pd.binu)<-colnames(pd.binu)<-colnames(comr.bin[[u]])
                        as.matrix(iCAMP::bmpd(comr.bin[[u]], pd.binu, abundance.weighted = weighted,unit.sum=unit.sum))
                      })
    bMPD.rand
  }
  
  # calculate across all samples #
  message("Now calculating observed betaMPD. Begin at ", date(),". Please wait...")
  gc()
  bMPD.obs<-lapply(1:bin.num,
                   function(i)
                   {
                     pd.bini=pd[pdid.bin[[i]],pdid.bin[[i]]]
                     gc()
                     rownames(pd.bini)<-colnames(pd.bini)<-colnames(com.bin[[i]])
                     as.matrix(iCAMP::bmpd(com.bin[[i]], pd.bini, abundance.weighted = weighted,unit.sum=unit.sum))
                   }) # calculate observed betaMPD.
  bMPD.obs<-array(unlist(bMPD.obs), dim = c(nrow(bMPD.obs[[1]]),ncol(bMPD.obs[[2]]),bin.num))
  gc()
  c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
  if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
  if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
  message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
  bMPD.rand<-parallel::parLapply(c1,1:rand,bMPD.random,
                                 pd.desc, pd.spname, pd.wd, pdid.bin,
                                 comms,sp.bin,bin.num,bin.lev,weighted,unit.sum,
                                 meta.group,meta.spool,meta.spall)
  parallel::stopCluster(c1)
  gc()
  bMPD.rand<-array(unlist(bMPD.rand),dim=c(nrow(bMPD.rand[[1]][[1]]),ncol(bMPD.rand[[1]][[1]]),length(bMPD.rand[[1]]),length(bMPD.rand)))
  
  if(sig.index=="RC")
  {
    bMPD.obsar=array(bMPD.obs,dim=dim(bMPD.rand))
    alpha1=apply(bMPD.obsar==bMPD.rand,c(1,2,3),sum)
    alpha=apply(bMPD.obsar>bMPD.rand,c(1,2,3),sum)
    alpha=(alpha+0.5*alpha1)/rand
    result=2*alpha-1
    gc()
    res.n=lapply(1:(dim(result)[3]), function(i){out=result[,,i];colnames(out)<-rownames(out)<-samp.name;out})
  }else if(sig.index %in% c("SES","bNRI")){
    bNRI=(bMPD.obs-apply(bMPD.rand,c(1,2,3),mean))/(apply(bMPD.rand,c(1,2,3),stats::sd))
    gc()
    res.n<-list()
    for(i in 1:(dim(bNRI)[3]))
    {
      diag(bNRI[,,i])<-0
      res.n[[i]]=bNRI[,,i]
      colnames(res.n[[i]])<-rownames(res.n[[i]])<-samp.name
    }
  }else if(sig.index=="Confidence"){
    bMPD.obsar=array(bMPD.obs,dim=dim(bMPD.rand))
    alpha=(apply(bMPD.obsar>bMPD.rand,c(1,2,3),sum))/rand
    alpha2=(apply(bMPD.obsar<bMPD.rand,c(1,2,3),sum))/rand
    alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
    result=alpha
    res.n=lapply(1:(dim(result)[3]), function(i){out=result[,,i];colnames(out)<-rownames(out)<-samp.name;out})
  }
  
  if(correct.special)
  {
    message("Now fixing special cases. Begin at ",date(),". Please wait...")
    sdm=(apply(bMPD.rand,c(1,2,3),stats::sd))
    correct.spec<-function(commi,resulti,pdi,sdmi,
                           meta.spool.bini,meta.freq.bini,meta.ab.bini,...)
    {
      # some special cases
      diag(sdmi)<-NA
      if(detail.null){special.ses=resulti;special.ses[]=0;special.rc<-special.conf<-special.ses}
      if(sum(sdmi==0,na.rm = TRUE)>0)
      {
        rownames(sdmi)<-colnames(sdmi)<-rownames(commi)
        sdc=dist.3col(sdmi)
        sdc0=sdc[which(sdc[,3]==0),,drop=FALSE]
        samp.sd0=unique(as.vector(as.matrix(sdc0[,1:2])))
        samp0=rownames(commi)[which(rowSums(commi)==0)]
        samp.ck=setdiff(samp.sd0,samp0)
        if(length(samp.ck)>0)
        {
          comm.ck=commi[which(rownames(commi) %in% samp.ck),,drop=FALSE]
          meta.group.ck=meta.group[which(rownames(meta.group) %in% samp.ck),,drop=FALSE]
          meta.spool.bini.ck=meta.spool.bini[match(unique(meta.group.ck[,1]),names(meta.spool))]
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
            nti.ck=iCAMP::NTI.cm(comm = comm.ck, dis = pdi, nworker = nworker,
                                 meta.group = meta.group.ck, meta.spool = meta.spool.bini.ck,
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
          
          if(special.method[1] %in% c("MPD", "NRI","both"))
          {
            nri.ck=iCAMP::NRI.cm(comm = comm.ck, dis = pdi, nworker = nworker,
                                 meta.group=meta.group.ck, meta.spool=meta.spool.bini.ck,
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
        
        samp1.id=which(rowSums(commi>0)==1)
        
        if(length(samp1.id)>0 | length(samp.sg)>0 | samp.sg.sum>0)
        {
          rcm=(iCAMP::RC.cm(comm=commi,rand=rand,na.zero=TRUE,nworker=nworker,
                            meta.group=meta.group,meta.frequency=meta.freq.bini,
                            meta.ab = meta.ab.bini,memory.G=memo.size.GB,weighted=weighted,
                            unit.sum=unit.sum,sig.index = 'RC',
                            silent=TRUE,dirichlet = dirichlet))$index
          
          if(length(samp1.id)>0)
          {
            requireNamespace("vegan")
            BCm=as.matrix(vegan::vegdist(commi,method = "euclidean",binary = TRUE))
            diag(BCm)=NA
            id.samesp=which(BCm==0,arr.ind = TRUE)
            id.same1=id.samesp[which((id.samesp[,1] %in% samp1.id)&(id.samesp[,2] %in% samp1.id)),,drop=FALSE]
            
            if(sig.index %in% c("RC","Confidence"))
            {
              resulti[id.same1]=(1-(2*(rcm[id.same1]<=0)))*1.1
            }else if(sig.index %in% c("SES","bNTI")){
              resulti[id.same1]=(1-(2*(rcm[id.same1]<=0)))*99
            }
            
            if(detail.null)
            {
              special.ses[id.same1]=(1-(2*(rcm[id.same1]<=0)))*99
              special.rc[id.same1]<-special.conf[id.same1]<-((1-(2*(rcm[id.same1]<=0)))*1.1)
            }
          }
          sd0.rcn=which(sdmi==0,arr.ind = TRUE)
          if(length(samp.sg)>0)
          {
            sd0.sg=sd0.rcn[which((rownames(sdmi)[sd0.rcn[,1]] %in% samp.sg)|(rownames(sdmi)[sd0.rcn[,2]] %in% samp.sg)),,drop=FALSE]
            if(sig.index %in% c("RC","Confidence"))
            {
              resulti[sd0.sg]=(1-(2*(rcm[sd0.sg]<=0)))*1.1
            }else if(sig.index %in% c("SES","bNTI")){
              resulti[sd0.sg]=(1-(2*(rcm[sd0.sg]<=0)))*99
            }
          }
          if(detail.null & samp.sg.sum>0)
          {
            sd0.sg.ses=sd0.rcn[which((rownames(sdmi)[sd0.rcn[,1]] %in% samp.sg.ses)|(rownames(sdmi)[sd0.rcn[,2]] %in% samp.sg.ses)),,drop=FALSE]
            special.ses[sd0.sg.ses]=(1-(2*(rcm[sd0.sg.ses]<=0)))*99
            sd0.sg.rc=sd0.rcn[which((rownames(sdmi)[sd0.rcn[,1]] %in% samp.sg.rc)|(rownames(sdmi)[sd0.rcn[,2]] %in% samp.sg.rc)),,drop=FALSE]
            special.rc[sd0.sg.rc]=(1-(2*(rcm[sd0.sg.rc]<=0)))*1.1
            sd0.sg.conf=sd0.rcn[which((rownames(sdmi)[sd0.rcn[,1]] %in% samp.sg.conf)|(rownames(sdmi)[sd0.rcn[,2]] %in% samp.sg.conf)),,drop=FALSE]
            special.conf[sd0.sg.conf]=(1-(2*(rcm[sd0.sg.conf]<=0)))*1.1
          }
        }
      }
      resulti[is.na(resulti)]=0
      if(detail.null){outii=list(resulti=resulti,special.ses=special.ses,special.rc=special.rc,special.conf=special.conf)}else{outii=resulti}
      outii
    }
    
    res.nn=lapply(1:length(res.n),
                 function(i)
                 {
                   sp.bini=rownames(sp.bin)[which(sp.bin[,1]==bin.lev[i])]
                   pdi=pd[pdid.bin[[i]],pdid.bin[[i]]]
                   rownames(pdi)<-colnames(pdi)<-colnames(com.bin[[i]])
                   meta.spool.bini=lapply(meta.spool,function(spj){spj[which(spj %in% sp.bini)]})
                   names(meta.spool.bini)=names(meta.spool)
                   meta.freq.bini=meta.frequency[,which(colnames(meta.frequency) %in% sp.bini),drop=FALSE]
                   meta.ab.bini=meta.ab[,which(colnames(meta.ab) %in% sp.bini),drop=FALSE]
                   correct.spec(commi = com.bin[[i]],resulti = res.n[[i]],
                                pdi = pdi, sdmi=sdm[,,i],meta.spool.bini=meta.spool.bini,
                                meta.freq.bini=meta.freq.bini,meta.ab.bini=meta.ab.bini)
                 })
    if(detail.null)
    {
      res.n=lapply(1:length(res.nn),function(i){res.nn[[i]]$resulti})
      spec.name2=iCAMP::dist.3col(res.nn[[1]]$special.ses)[,1:2,drop=FALSE]
      xxspec<-function(ind,...)
      {
        outx=sapply(1:length(res.nn),function(i){iCAMP::dist.3col(res.nn[[i]][[ind]])[,3]})
        colnames(outx)=paste0("bin",1:ncol(outx))
        outx
      }
      special.ses.bin=data.frame(spec.name2,xxspec("special.ses"),stringsAsFactors = FALSE)
      special.rc.bin=data.frame(spec.name2,xxspec("special.rc"),stringsAsFactors = FALSE)
      special.conf.bin=data.frame(spec.name2,xxspec("special.conf"),stringsAsFactors = FALSE)
    }else{
      res.n=res.nn
    }
  }
  
  output=list(index=res.n,sig.index=sig.index)
  if(output.bMPD)
  {
    colnames(bMPD.obs)<-rownames(bMPD.obs)<-samp.name
    output=c(output,list(betaMPD.obs=lapply(1:(dim(bMPD.obs)[3]),function(i){bMPD.obs[,,i]})))
  }
  if(detail.null)
  {
    rownames(bMPD.rand)=samp.name
    colnames(bMPD.rand)=samp.name
    samp2.name=(iCAMP::dist.3col(bMPD.rand[,,1,1]))[,1:2,drop=FALSE]
    bMPD.randm=lapply(1:(dim(bMPD.rand)[3]),
                      function(i)
                      {
                        outi=matrix(sapply(1:(dim(bMPD.rand)[4]),
                                    function(j)
                                    {
                                      (iCAMP::dist.3col(bMPD.rand[,,i,j]))[,3]
                                    }),ncol=(dim(bMPD.rand)[4]))
                        colnames(outi)=paste0("rand",1:ncol(outi))
                        data.frame(samp2.name,outi,stringsAsFactors = FALSE)
                      })
    if(correct.special)
    {
      special.crct=list(special.ses.bin=special.ses.bin,special.rc.bin=special.rc.bin,special.conf.bin=special.conf.bin)
    }else{
      special.crct=NULL
    }
    output=c(output,list(rand=bMPD.randm, special.crct=special.crct))
  }
  output
}