RC.bin.cm<-function(com,sp.bin,rand=1000,na.zero=TRUE,
                    meta.group=NULL,meta.frequency=NULL,meta.ab=NULL,
                    nworker=4,memory.G=50,big.method=c("loop","no"),
                    weighted=TRUE,unit.sum=NULL,
                    sig.index=c("RC","Confidence","SES"),
                    detail.null=FALSE,output.bray=FALSE,
                    taxo.metric="bray", transform.method=NULL,
                    logbase=2, dirichlet=FALSE)
{
  requireNamespace("vegan")
  requireNamespace("parallel")
  
  if(max(rowSums(com,na.rm = TRUE))<=1 & (!dirichlet))
  {
    warning("The values in com are less than 1, thus considered as proportional data, Dirichlet distribution is used to assign abundance in null model.")
    dirichlet=TRUE
  }
  
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memory.G*1024)
    {
      memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
      if(inherits(memotry,"try-error")){warning(memotry[1])}
    }
  }
  
  if(!weighted){com[com>0]=1}
  
  ##########################
  ## check IDs in meta communities
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
  if(is.null(unit.sum))
  {
    com.ra=comm/rowSums(comm)
    com.ra[which(rowSums(comm)==0),]=0
  }else{
    com.ra=comm/unit.sum
    com.ra[which(unit.sum==0),]=0
  }
  comras=lapply(meta.lev,function(mi){sampi=rownames(meta.group)[which(meta.group[,1]==mi)];comrai=com.ra[which(rownames(com.ra) %in% sampi),,drop=FALSE];comrai[,colSums(comrai)>0,drop=FALSE]})
  
  if(!is.null(meta.frequency))
  {
    if(sum(!(colnames(comm) %in% colnames(meta.frequency)))){stop('comm has some species not included in meta.frequence.')}
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
  
  if(!is.null(meta.ab))
  {
    if(sum(!(colnames(comm) %in% colnames(meta.ab)))){stop('comm has some species not included in meta.ab')}
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
      meta.ab[i,match(colnames(comms[[i]]),colnames(comm))]=colMeans(comras[[i]])
    }
  }
  
  ########################
  sig.index=sig.index[1]
  if(!(sig.index %in% c("RC","Confidence","SES"))){stop("wrong sig.index for RC.pc.")}
  
  bin.lev=levels(as.factor(sp.bin[,1]))
  bin.num=length(bin.lev)
  
  
  bin.BC<-function(comi,bin.lev,bin.num,sp.bin,
                   unit.sum,na.zero,
                   taxo.metric,transform.method,logbase)
  {
    if(!is.null(transform.method))
    {
      if(inherits(transform.method,"function"))
      {
        comit=transform.method(comi)
      }else{
        comit=vegan::decostand(comi,method = transform.method,logbase = logbase,na.rm=TRUE)
      }
    }else{comit=comi}
    
    com.bin=lapply(1:bin.num, function(i){comit[,match(rownames(sp.bin)[which(sp.bin==bin.lev[i])],colnames(comit))]})
    
    lapply(1:bin.num,
           function(i)
           {
             combini=com.bin[[i]]
             if(taxo.metric=='bray' & is.null(transform.method))
             {
               if(is.null(unit.sum))
               {
                 combinit=combini/rowSums(combini)
                 combinit[which(rowSums(combini)==0),]=0
                 binBC=vegan::vegdist(combini,method="bray")
               }else{
                 combinit=combini/unit.sum
                 combinit[which(unit.sum==0),]=0
                 binBC=vegan::vegdist(combinit,method="manhattan")/2
               }
             }else{
               binBC=vegan::vegdist(combini,method = taxo.metric)
             }
             if(na.zero){binBC[is.na(binBC)]=0}
             as.matrix(binBC)
           })
  }
  
  BC.obs<-bin.BC(comi=com,bin.lev=bin.lev,bin.num=bin.num,sp.bin=sp.bin,
                 unit.sum=unit.sum,na.zero=na.zero,taxo.metric=taxo.metric,
                 transform.method=transform.method,logbase=logbase)
  BC.obs<-array(unlist(BC.obs),dim = c(nrow(BC.obs[[1]]),ncol(BC.obs[[1]]),bin.num))
  
  #######################
  Si=lapply(comms,function(comi){rowSums(comi>0)})
  Ni=lapply(comms,function(comi){rowSums(comi)})
  
  BC.rand<-function(j,meta.group,meta.frequency,meta.ab,Si,Ni,na.zero,unit.sum,
                    taxo.metric,transform.method,logbase,dirichlet,
                    sp.bin,bin.num,bin.BC,bin.lev)
  {
    requireNamespace("vegan")
    com.rd=matrix(0,nrow=nrow(meta.group),ncol=ncol(meta.frequency))
    rownames(com.rd)=rownames(meta.group)
    colnames(com.rd)=colnames(meta.frequency)
    for(k in 1:length(Si))
    {
      meta.freqk=meta.frequency[k,]
      meta.abk=meta.ab[k,]
      Sik=Si[[k]]
      Nik=Ni[[k]]
      if(!dirichlet)
      {
        for(i in 1:length(Sik))
        {
          if(Sik[i]!=0)
          {
            id.sp<-sample(1:ncol(com.rd),Sik[i],replace=FALSE,prob=meta.freqk)
            if(length(id.sp)==1){count=rep(id.sp,Nik[i]-1)}else{
              count<-sample(id.sp,(Nik[i]-Sik[i]),replace=TRUE,prob=meta.abk[id.sp])
            }
            tabk<-table(count)
            idik=match(names(Sik)[i],rownames(com.rd))
            com.rd[idik,as.numeric(names(tabk))]=as.vector(tabk)
            com.rd[idik,id.sp]=com.rd[idik,id.sp]+1
          }
        }
      }else{
        for(i in 1:length(Sik))
        {
          if(Sik[i]!=0)
          {
            id.sp<-sample(1:ncol(com.rd),Sik[i],replace=FALSE,prob=meta.freqk)
            idik=match(names(Sik)[i],rownames(com.rd))
            if(length(id.sp)==1){com.rd[idik,id.sp]=1}else{
              requireNamespace("DirichletReg")
              com.rd[idik,id.sp]=DirichletReg::rdirichlet(n=1,alpha = meta.abk[id.sp])
            }
          }
        }
      }
    }
    
    BCrand<-bin.BC(comi=com.rd,bin.lev=bin.lev,bin.num=bin.num,sp.bin=sp.bin,
                   unit.sum=unit.sum,na.zero=na.zero,taxo.metric=taxo.metric,
                   transform.method=transform.method,logbase=logbase)
    BCrand
  }
  
  if(big.method[1]=="loop")
  {
    if(detail.null){BC.rdl=list()}
    if(sig.index=="RC")
    {
      alpha=array(0,dim=dim(BC.obs))
    }else if(sig.index=="Confidence"){
      alpha<-alpha2<-array(0,dim=dim(BC.obs))
    }
      
    for(j in 1:rand)
    {
      message("Now rand BC j=",j," in ",rand,". ",date())
      BC.rd=BC.rand(j,meta.group,meta.frequency,meta.ab,Si,Ni,na.zero,unit.sum,
                    taxo.metric,transform.method,logbase,dirichlet,
                    sp.bin,bin.num,bin.BC,bin.lev)
      BCrd.array=array(unlist(BC.rd),dim=c(nrow(BC.rd[[1]]),ncol(BC.rd[[1]]),length(BC.rd)))
      if(detail.null | (sig.index=="SES")){BC.rdl[[j]]=BCrd.array}
      if(sig.index=="RC")
      {
        idxx=which(BCrd.array<BC.obs)
        alpha[idxx]=alpha[idxx]+1
        idxx=which(BCrd.array==BC.obs)
        alpha[idxx]=alpha[idxx]+0.5
      }else if(sig.index=="Confidence"){
        idxx=which(BC.obs>BCrd.array,arr.ind = TRUE)
        alpha[idxx]=alpha[idxx]+1
        idxx2=which(BC.obs<BCrd.array,arr.ind = TRUE)
        alpha2[idxx]=alpha2[idxx]+1
      }
      gc()
    }
    message("----now calculating sig.index of Bray at ",date(),"----")
    if(sig.index=="RC")
    {
      alpha=alpha/rand
      rc=(alpha-0.5)*2
      result=list()
      for(i in 1:(dim(rc)[3]))
      {
        result[[i]]=rc[,,i]
        colnames(result[[i]])<-rownames(result[[i]])<-rownames(com)
      }
    }else if(sig.index=="Confidence"){
      alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
      colnames(alpha)<-rownames(alpha)<-rownames(com)
      result=lapply(1:(dim(alpha)[3]),function(i){alpha[,,i]})
    }else if(sig.index=="SES"){
      result=lapply(1:(dim(BC.obs)[3]),
                    function(i)
                    {
                      BC.rdli=sapply(1:length(BC.rdl),function(x){BC.rdl[[x]][,,i]},simplify = "array")
                      SESi=(BC.obs[,,i]-apply(BC.rdli,c(1,2),mean))/(apply(BC.rdli,c(1,2),stats::sd))
                      diag(SESi)<-0
                      colnames(SESi)<-rownames(SESi)<-rownames(com)
                      SESi
                    })
    }
  }else{
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
    message("Now parallel computing. begin at ", date(),". Please wait...")
    BC.rd<-parallel::parLapply(c1,1:rand,BC.rand,
                               meta.group,meta.frequency,meta.ab,Si,Ni,na.zero,unit.sum,
                               taxo.metric,transform.method,logbase,dirichlet,
                               sp.bin,bin.num,bin.BC,bin.lev)
    parallel::stopCluster(c1)
    BC.rd=array(unlist(BC.rd),dim=c(nrow(BC.rd[[1]][[1]]),ncol(BC.rd[[1]][[1]]),length(BC.rd[[1]]),length(BC.rd)))
    gc()
    message("----now calculating sig.index of Bray at ",date(),"----")
    if(sig.index=="RC")
    {
      comp<-function(x,c){(x<c)+0.5*(x==c)}
      alpha=array(rowSums(apply(BC.rd,4,comp,c=BC.obs)),dim=dim(BC.obs))/rand
      rc=(alpha-0.5)*2
      result=list()
      for(i in 1:(dim(rc)[3]))
      {
        result[[i]]=rc[,,i]
        colnames(result[[i]])<-rownames(result[[i]])<-rownames(com)
      }
    }else if(sig.index=="Confidence"){
      BC.obsar=array(BC.obs,dim=dim(BC.rd))
      alpha=(apply(BC.obsar>BC.rd,c(1,2,3),sum))/rand
      alpha2=(apply(BC.obsar<BC.rd,c(1,2,3),sum))/rand
      alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
      result=lapply(1:(dim(alpha)[3]), function(i){out=alpha[,,i];colnames(out)<-rownames(out)<-rownames(com);out})
    }else if(sig.index=="SES"){
      SES=(BC.obs-apply(BC.rd,c(1,2,3),mean))/(apply(BC.rd,c(1,2,3),stats::sd))
      result=lapply(1:(dim(SES)[3]), function(i){out=SES[,,i];diag(out)=0;colnames(out)<-rownames(out)<-rownames(com);out})
    }
  }
  
  output=list(index=result,sig.index=sig.index)
  if(output.bray)
  {
    colnames(BC.obs)<-rownames(BC.obs)<-rownames(com)
    output=c(output,list(BC.obs=lapply(1:(dim(BC.obs)[3]),function(i){BC.obs[,,i]})))
  }
  if(detail.null)
  {
    if(big.method[1]=="loop")
    {
      for(i in 1:length(BC.rdl)){colnames(BC.rdl[[i]])<-rownames(BC.rdl[[i]])<-rownames(com)}
      samp2.name=(iCAMP::dist.3col(BC.rdl[[1]][,,1]))[,1:2,drop=FALSE]
      BC.randm=lapply(1:(dim(BC.rdl[[1]])[3]),
                      function(i)
                      {
                        outi=matrix(sapply(1:length(BC.rdl),
                                    function(j)
                                    {
                                      (iCAMP::dist.3col(BC.rdl[[j]][,,i]))[,3]
                                    }),ncol=length(BC.rdl))
                        colnames(outi)=paste0("rand",1:ncol(outi))
                        data.frame(samp2.name,outi,stringsAsFactors = FALSE)
                      })
    }else{
      rownames(BC.rd)<-colnames(BC.rd)<-rownames(com)
      samp2.name=(iCAMP::dist.3col(BC.rd[,,1,1]))[,1:2,drop=FALSE]
      BC.randm=lapply(1:(dim(BC.rd)[3]),
                         function(i)
                         {
                           outi=matrix(sapply(1:(dim(BC.rd)[4]),
                                       function(j)
                                       {
                                         (iCAMP::dist.3col(BC.rd[,,i,j]))[,3]
                                       }),ncol=(dim(BC.rd)[4]))
                           colnames(outi)=paste0("rand",1:ncol(outi))
                           data.frame(samp2.name,outi,stringsAsFactors = FALSE)
                         })
    }
    output=c(output,list(rand=BC.randm))
  }
  output
}