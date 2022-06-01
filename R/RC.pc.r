RC.pc<-function(comm,rand=1000,na.zero=TRUE,nworker=4,
                memory.G=50,weighted=TRUE,unit.sum=NULL,
                meta.ab=NULL,sig.index=c("RC","Confidence","SES"),
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
  com<-comm[,colSums(comm)>0]
  if(!weighted){com[com>0]=1}
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
  
  BC.obs<-BCnew(com,unit.sum,na.zero,
                taxo.metric,transform.method,logbase)
  
  if(is.null(meta.ab))
  {
    if(is.null(unit.sum))
    {
      com.ra=com/rowSums(com)
      com.ra[which(rowSums(com)==0),]=0
    }else{
      com.ra=com/unit.sum
      com.ra[which(unit.sum==0),]=0
    }
    prob.ab=colMeans(com.ra)
  }else{
    if(is.null(names(meta.ab)))
    {
      if(length(meta.ab)!=ncol(comm)){stop("meta.ab setting is wrong.")}else{prob.ab=meta.ab}
    }else{
      prob.ab=rep(0,ncol(comm))
      meta.ab=meta.ab[which(names(meta.ab) %in% colnames(comm))]
      prob.ab[match(names(meta.ab),colnames(comm))]=meta.ab
    }
  }
  
  com.rd0=com
  com.rd0[]=0
  id<-(1:ncol(com))
  prob.sp<-colSums(com>0)
  #prob.ab<-colSums(com)
  Si<-rowSums(com>0)
  Ni<-rowSums(com)
  samp.num=nrow(com)
  
  BC.rand<-function(j,com.rd0,samp.num,id,prob.sp,prob.ab,Si,Ni,na.zero,unit.sum,
                    BCnew,taxo.metric,transform.method,logbase,dirichlet)
  {
    requireNamespace("vegan")
    
    com.rd=com.rd0
    if(!dirichlet)
    {
      for(i in 1:samp.num)
      {
        if(Si[i]==0){com.rd[i,]=0}else{
          id.sp<-sample(id,Si[i],replace=FALSE,prob=prob.sp)
          if(length(id.sp)==1){count=rep(id.sp,Ni[i])}else{
            count<-sample(id.sp,(Ni[i]-Si[i]),replace=TRUE,prob=prob.ab[id.sp])
          }
          table<-table(count)
          com.rd[i,as.numeric(names(table))]=as.vector(table)
          com.rd[i,id.sp]=com.rd[i,id.sp]+1
        }
      }
    }else{
      for(i in 1:samp.num)
      {
        if(Si[i]==0){com.rd[i,]=0}else{
          id.sp<-sample(id,Si[i],replace=FALSE,prob=prob.sp)
          if(length(id.sp)==1){com.rd[i,id.sp]=1}else{
            requireNamespace("DirichletReg")
            com.rd[i,id.sp]=DirichletReg::rdirichlet(n=1,alpha = prob.ab[id.sp])
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
  BC.rd<-parallel::parLapply(c1,1:rand,BC.rand,com.rd0=com.rd0,samp.num=samp.num,id=id,
                             prob.sp=prob.sp,prob.ab=prob.ab,Si=Si,Ni=Ni,na.zero=na.zero,unit.sum=unit.sum,
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