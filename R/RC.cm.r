RC.cm<-function(comm,rand=1000,na.zero=TRUE,nworker=4,
                meta.group=NULL,meta.frequency=NULL,meta.ab=NULL,
                memory.G=50,weighted=TRUE,unit.sum=NULL,
                sig.index=c("RC","Confidence","SES"),
                detail.null=FALSE,output.bray=FALSE,silent=FALSE,
                taxo.metric="bray", transform.method=NULL, logbase=2,
                dirichlet=FALSE)
{
  # v20200728 add sig.index, detail.null, output.bray
  
  requireNamespace("vegan")
  requireNamespace("parallel")
  if(max(rowSums(comm,na.rm = TRUE))<=1 & (!dirichlet))
  {
    warning("The values in comm are less than 1, thus considered as proportional data, Dirichlet distribution is used to assign abundance in null model.")
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
  if(!weighted){comm[comm>0]=1}
  ####################
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
  
  ####################
  ## calculate observed beta
  sig.index=sig.index[1]
  if(!(sig.index %in% c("RC","Confidence","SES"))){stop("wrong sig.index for RC.pc.")}
  if(!is.null(unit.sum)){if(length(unit.sum)==1){unit.sum=rep(unit.sum,nrow(comm))}}
  
  BCnew<-function(comi,unit.sum,na.zero,
                  taxo.metric,transform.method,logbase)
  {
    if(taxo.metric=='bray' & is.null(transform.method))
    {
      if(is.null(unit.sum))
      {
        comit=comi/rowSums(comi)
        comit[which(rowSums(comi)==0),]=0
        BC=vegan::vegdist(comit,method="bray")
      }else{
        comit=comi/unit.sum
        comit[which(unit.sum==0),]=0
        BC=vegan::vegdist(comit,method="manhattan")/2
      }
    }else{
      if(!is.null(transform.method))
      {
        if(inherits(transform.method,"function"))
        {
          comit=transform.method(comi)
        }else{
          comit=vegan::decostand(comi,method = transform.method,logbase = logbase,na.rm=TRUE)
        }
      }else{comit=comi}
      BC=vegan::vegdist(comit,method = taxo.metric)
    }
    if(na.zero){BC[is.na(BC)]=0}
    as.matrix(BC)
  }
  
  BC.obs<-BCnew(comm,unit.sum,na.zero,
                taxo.metric,transform.method,logbase)
  
  #########################
  # Null model
  
  Si=lapply(comms,function(comi){rowSums(comi>0)})
  Ni=lapply(comms,function(comi){rowSums(comi)})
  
  
  BC.rand<-function(j,meta.group,meta.frequency,meta.ab,Si,Ni,na.zero,unit.sum,
                    BCnew,taxo.metric,transform.method,logbase,dirichlet)
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
    
    BCnew(com.rd,unit.sum,na.zero,
          taxo.metric,transform.method,logbase)
  }
  
  c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
  if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
  if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
  if(!silent){message("Now parallel computing. begin at ", date(),". Please wait...")}
  BC.rd<-parallel::parLapply(c1,1:rand,BC.rand,meta.group=meta.group,meta.frequency=meta.frequency,meta.ab=meta.ab,
                             Si=Si,Ni=Ni,na.zero=na.zero,unit.sum=unit.sum,
                             BCnew=BCnew,taxo.metric=taxo.metric,
                             transform.method=transform.method,logbase=logbase,dirichlet=dirichlet)
  parallel::stopCluster(c1)
  
  BC.rd=array(unlist(BC.rd),dim=c(nrow(BC.rd[[1]]),ncol(BC.rd[[1]]),length(BC.rd)))
  gc()
  
  if(!silent){message("----now calculating sig.index of Bray at ",date(),"----")}
  if(sig.index=="RC")
  {
    comp<-function(x,c){(x<c)+0.5*(x==c)}
    alpha=matrix(rowSums(apply(BC.rd,3,comp,c=BC.obs)),nrow=nrow(BC.obs))/rand
    result=(alpha-0.5)*2
  }else if(sig.index=="SES"){
    SES=(BC.obs-apply(BC.rd,c(1,2),mean))/(apply(BC.rd,c(1,2),stats::sd))
    diag(SES)<-0
    result=SES
  }else if(sig.index=="Confidence"){
    BC.obsar=array(BC.obs,dim=dim(BC.rd))
    alpha=(apply(BC.obsar>BC.rd,c(1,2),sum))/rand
    alpha2=(apply(BC.obsar<BC.rd,c(1,2),sum))/rand
    alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
    result=alpha
  }
  
  rownames(result)=rownames(BC.obs)
  colnames(result)=colnames(BC.obs)
  
  output=list(index=result,sig.index=sig.index)
  if(output.bray)
  {
    output=c(output,list(BC.obs=BC.obs))
  }
  if(detail.null)
  {
    rownames(BC.rd)=rownames(BC.obs)
    colnames(BC.rd)=colnames(BC.obs)
    BC.randm=matrix(sapply(1:(dim(BC.rd)[3]),
                       function(i)
                       {
                         (iCAMP::dist.3col(BC.rd[,,i]))[,3]
                       }),ncol=(dim(BC.rd)[3]))
    colnames(BC.randm)=paste0("rand",1:ncol(BC.randm))
    output=c(output,list(rand=data.frame((iCAMP::dist.3col(BC.rd[,,1]))[,1:2,drop=FALSE],
                                         BC.randm,stringsAsFactors = FALSE)))
  }
  output
}
