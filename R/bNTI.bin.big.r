bNTI.bin.big<-function(comm, pd.desc, pd.spname, pd.wd, pdid.bin, sp.bin,
                       spname.check=FALSE,nworker=4, memo.size.GB=50,
                       weighted=c(TRUE,FALSE),rand=1000,output.bMNTD=c(FALSE,TRUE),
                       sig.index=c("SES","Confidence","RC","bNTI"),unit.sum=NULL,
                       correct.special=FALSE,detail.null=FALSE,
                       special.method=c("MNTD","MPD","both"),
                       ses.cut=1.96, rc.cut=0.95, conf.cut=0.975,
                       exclude.conspecifics = FALSE)
{
  #v20200728 add conf.cut, detail.null. change RC to sig.index.
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
  weighted=weighted[1];output.bMNTD=output.bMNTD[1]
  sig.index=sig.index[1]
  special.method=special.method[1]
  if(!(sig.index %in% c("SES","Confidence","RC","bNTI"))){stop("wrong sig.index for bNTI.bin.big.")}
  
  requireNamespace("bigmemory")
  pd=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
  
  # match
  if(spname.check)
  {
    sp.check=iCAMP::match.name(name.check=pd.spname,cn.list = list(comm=comm))
    comm=sp.check$comm
  }
  samp.name=rownames(comm)
  sp.name=colnames(comm)
  sp.num=ncol(comm)
  
  # bin comm and dis
  bin.lev=levels(as.factor(sp.bin[,1]))
  bin.num=length(bin.lev)
  com.bin=lapply(1:bin.num, function(i){comm[,match(rownames(sp.bin)[which(sp.bin==i)],colnames(comm))]})
  #dis.bin=lapply(1:bin.num, function(i){id=match(rownames(sp.bin)[which(sp.bin==i)],rownames(dis));dis[id,id]})
  
  ## randomization function ##
  requireNamespace("permute")
  perm=permute::shuffleSet(sp.num,rand)
  if(sp.num==2){perm=matrix(c(2,1),nrow=1)}
  if(nrow(perm)<rand)
  {
    perm=rbind(perm,1:sp.num)
    rand=nrow(perm)
  }
  
  bMNTD.random<-function(i,pd.desc, pd.spname, pd.wd, pdid.bin,
                         com,sp.bin,bin.num,weighted,perm,
                         unit.sum,exclude.conspecifics)
  {
    requireNamespace("iCAMP")
    com.rand=com[,perm[i,]]
    colnames(com.rand)<-colnames(com)
    comr.bin=lapply(1:bin.num, function(j){com.rand[,match(rownames(sp.bin)[which(sp.bin==j)],colnames(com.rand))]})
    gc()
    requireNamespace("bigmemory")
    pd=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
    bMNTD.rand<-lapply(1:bin.num,
                       function(u)
                       {
                         pd.binu=pd[pdid.bin[[u]],pdid.bin[[u]]]
                         gc()
                         rownames(pd.binu)<-colnames(pd.binu)<-colnames(comr.bin[[u]])
                         as.matrix(iCAMP::bmntd(comr.bin[[u]], pd.binu,
                                                abundance.weighted = weighted,unit.sum=unit.sum,
                                                exclude.conspecifics=exclude.conspecifics))
                       })
    bMNTD.rand
  }
  
  
  # calculate across all samples #
  message("Now calculating observed betaMNTD. Begin at ", date(),". Please wait...")
  gc()
  bMNTD.obs<-lapply(1:bin.num,
                    function(i)
                    {
                      pd.bini=pd[pdid.bin[[i]],pdid.bin[[i]]]
                      gc()
                      rownames(pd.bini)<-colnames(pd.bini)<-colnames(com.bin[[i]])
                      as.matrix(iCAMP::bmntd(com.bin[[i]], pd.bini, abundance.weighted = weighted,
                                             unit.sum=unit.sum,exclude.conspecifics=exclude.conspecifics))
                    }) # calculate observed betaMNTD.
  bMNTD.obs<-array(unlist(bMNTD.obs), dim = c(nrow(bMNTD.obs[[1]]),ncol(bMNTD.obs[[1]]),bin.num))
  gc()
  c1<-parallel::makeCluster(nworker,type="PSOCK")
  message("Now randomizing by parallel computing. Begin at ", date(),". Please wait...")
  bMNTD.rand<-parallel::parLapply(c1,1:rand,bMNTD.random,pd.desc=pd.desc, pd.spname=pd.spname,
                        pd.wd=pd.wd, pdid.bin=pdid.bin, com=comm,sp.bin=sp.bin,
                        bin.num=bin.num,weighted=weighted,perm=perm,
                        unit.sum=unit.sum,exclude.conspecifics=exclude.conspecifics)
  parallel::stopCluster(c1)
  gc()
  bMNTD.rand<-array(unlist(bMNTD.rand),dim=c(nrow(bMNTD.rand[[1]][[1]]),ncol(bMNTD.rand[[1]][[1]]),length(bMNTD.rand[[1]]),length(bMNTD.rand)))
  
  if(sig.index=="RC")
  {
    bMNTD.obsar=array(bMNTD.obs,dim=dim(bMNTD.rand))
    alpha1=apply(bMNTD.obsar==bMNTD.rand,c(1,2,3),sum)
    alpha=apply(bMNTD.obsar>bMNTD.rand,c(1,2,3),sum)
    alpha=(alpha+0.5*alpha1)/rand
    result=2*alpha-1
    gc()
    res.n=lapply(1:(dim(result)[3]), function(i){out=result[,,i];colnames(out)<-rownames(out)<-samp.name;out})
  }else if(sig.index %in% c("SES","bNRI")){
    bNTI=(bMNTD.obs-apply(bMNTD.rand,c(1,2,3),mean))/(apply(bMNTD.rand,c(1,2,3),stats::sd))
    gc()
    res.n<-list()
    for(i in 1:(dim(bNTI)[3]))
    {
      diag(bNTI[,,i])<-0
      res.n[[i]]=bNTI[,,i]
      colnames(res.n[[i]])<-rownames(res.n[[i]])<-samp.name
    }
  }else if(sig.index=="Confidence"){
    bMNTD.obsar=array(bMNTD.obs,dim=dim(bMNTD.rand))
    alpha=(apply(bMNTD.obsar>bMNTD.rand,c(1,2,3),sum))/rand
    alpha2=(apply(bMNTD.obsar<bMNTD.rand,c(1,2,3),sum))/rand
    alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
    result=alpha
    res.n=lapply(1:(dim(result)[3]), function(i){out=result[,,i];colnames(out)<-rownames(out)<-samp.name;out})
  }
  
  if(correct.special)
  {
    message("Now fixing special cases. Begin at ",date(),". Please wait...")
    sdm=(apply(bMNTD.rand,c(1,2,3),stats::sd))
    correct.spec<-function(commi,resulti,pdi,sdmi,...)
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
          samp.sg=character(0)
          if(detail.null)
          {
            samp.sg.ses<-samp.sg.rc<-samp.sg.conf<-character(0)
            sig.nti="all"
          }else{
            if(sig.index %in% c("bNTI","SES")){sig.nti="SES";nti.cut=ses.cut}else if(sig.index=="RC"){sig.nti="RC";nti.cut=rc.cut}else if(sig.index=="Confidence"){sig.nti="Confidence";nti.cut=conf.cut}
          }
          if(special.method[1] %in% c("MNTD","NTI","both"))
          {
            nti.ck=iCAMP::NTI.p(comm = comm.ck, dis = pdi, nworker = nworker,
                         memo.size.GB = memo.size.GB,weighted = weighted,
                         rand = rand,output.MNTD = FALSE,sig.index = sig.nti)
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
            nri.ck=iCAMP::NRI.p(comm = comm.ck, dis = pdi, nworker = nworker,
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
        }else{
          samp.sg=character(0)
          samp.sg.sum=0
        }
        
        samp1.id=which(rowSums(commi>0)==1)
        
        if(length(samp1.id)>0 | length(samp.sg)>0 | samp.sg.sum>0)
        {
          rcm=(iCAMP::RC.pc(comm=commi,rand=rand,na.zero=TRUE,nworker=nworker,
                            memory.G=memo.size.GB,weighted=weighted,
                            unit.sum=unit.sum,silent=TRUE))$index
          
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
                   pdi=pd[pdid.bin[[i]],pdid.bin[[i]]]
                   rownames(pdi)<-colnames(pdi)<-colnames(com.bin[[i]])
                   correct.spec(commi = com.bin[[i]],resulti = res.n[[i]],
                                pdi = pdi, sdmi=sdm[,,i])
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
  if(output.bMNTD)
  {
    colnames(bMNTD.obs)<-rownames(bMNTD.obs)<-samp.name
    output=c(output,list(betaMNTD.obs=lapply(1:(dim(bMNTD.obs)[3]),function(i){bMNTD.obs[,,i]})))
  }
  if(detail.null)
  {
    rownames(bMNTD.rand)=samp.name
    colnames(bMNTD.rand)=samp.name
    samp2.name=(iCAMP::dist.3col(bMNTD.rand[,,1,1]))[,1:2,drop=FALSE]
    bMNTD.randm=lapply(1:(dim(bMNTD.rand)[3]),
                      function(i)
                      {
                        outi=sapply(1:(dim(bMNTD.rand)[4]),
                                    function(j)
                                    {
                                      (iCAMP::dist.3col(bMNTD.rand[,,i,j]))[,3]
                                    })
                        colnames(outi)=paste0("rand",1:ncol(outi))
                        data.frame(samp2.name,outi,stringsAsFactors = FALSE)
                      })
    if(correct.special)
    {
      special.crct=list(special.ses.bin=special.ses.bin,special.rc.bin=special.rc.bin,special.conf.bin=special.conf.bin)
    }else{special.crct=NULL}
    output=c(output,list(rand=bMNTD.randm, special.crct=special.crct))
  }
  
  output
}