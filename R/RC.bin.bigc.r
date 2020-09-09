RC.bin.bigc<-function(com,sp.bin,rand=1000,na.zero=TRUE,
                      nworker=4,memory.G=50,big.method=c("loop","no"),
                      weighted=TRUE,unit.sum=NULL,meta.ab=NULL,
                      sig.index=c("RC","Confidence","SES"),
                      detail.null=FALSE,output.bray=FALSE)
{
  requireNamespace("vegan")
  requireNamespace("parallel")
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memory.G*1024)
    {
      memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
      if(class(memotry)=="try-error"){warning(memotry[1])}
    }
  }
  if(!weighted){com[com>0]=1}
  sig.index=sig.index[1]
  if(!(sig.index %in% c("RC","Confidence","SES"))){stop("wrong sig.index for RC.pc.")}
  
  bin.lev=levels(as.factor(sp.bin[,1]))
  bin.num=length(bin.lev)
  com.bin=lapply(1:bin.num, function(i){com[,match(rownames(sp.bin)[which(sp.bin==i)],colnames(com))]})
  
  bin.BC<-function(combini,unit.sum,na.zero)
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
    if(na.zero){binBC[is.na(binBC)]=0}
    as.matrix(binBC)
  }
  
  BC.obs<-lapply(1:bin.num,function(i){bin.BC(com.bin[[i]],unit.sum = unit.sum,na.zero = na.zero)})
  BC.obs<-array(unlist(BC.obs),dim = c(nrow(BC.obs[[1]]),ncol(BC.obs[[1]]),bin.num))
  
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
      if(length(meta.ab)!=ncol(com)){stop("meta.ab setting is wrong.")}else{prob.ab=meta.ab}
    }else{
      prob.ab=rep(0,ncol(com))
      meta.ab=meta.ab[which(names(meta.ab) %in% colnames(com))]
      prob.ab[match(names(meta.ab),colnames(com))]=meta.ab
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
  
  BC.rand<-function(j,com.rd0,samp.num,id,prob.sp,prob.ab,Si,Ni,na.zero,sp.bin,bin.num,unit.sum)
  {
    requireNamespace("vegan")
    bin.BC<-function(combini,unit.sum,na.zero)
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
      if(na.zero){binBC[is.na(binBC)]=0}
      as.matrix(binBC)
    }
    com.rd=com.rd0
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
    comr.bin=lapply(1:bin.num, function(m){com.rd[,match(rownames(sp.bin)[which(sp.bin==m)],colnames(com.rd))]})
    BCrand<-lapply(1:bin.num,function(m){bin.BC(comr.bin[[m]],unit.sum = unit.sum,na.zero = na.zero)})
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
      BC.rd=BC.rand(j,com.rd0=com.rd0,samp.num=samp.num,id=id,prob.sp=prob.sp,prob.ab=prob.ab,Si=Si,Ni=Ni,na.zero=na.zero,sp.bin=sp.bin,bin.num=bin.num,unit.sum=unit.sum)
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
    c1<-parallel::makeCluster(nworker,type="PSOCK")
    message("Now parallel computing. begin at ", date(),". Please wait...")
    BC.rd<-parallel::parLapply(c1,1:rand,BC.rand,com.rd0=com.rd0,samp.num=samp.num,
                     id=id,prob.sp=prob.sp,prob.ab=prob.ab,Si=Si,Ni=Ni,
                     na.zero=na.zero,sp.bin=sp.bin,bin.num=bin.num,unit.sum=unit.sum)
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
                        outi=sapply(1:length(BC.rdl),
                                    function(j)
                                    {
                                      (iCAMP::dist.3col(BC.rdl[[j]][,,i]))[,3]
                                    })
                        colnames(outi)=paste0("rand",1:ncol(outi))
                        data.frame(samp2.name,outi,stringsAsFactors = FALSE)
                      })
    }else{
      rownames(BC.rd)<-colnames(BC.rd)<-rownames(com)
      samp2.name=(iCAMP::dist.3col(BC.rd[,,1,1]))[,1:2,drop=FALSE]
      BC.randm=lapply(1:(dim(BC.rd)[3]),
                         function(i)
                         {
                           outi=sapply(1:(dim(BC.rd)[4]),
                                       function(j)
                                       {
                                         (iCAMP::dist.3col(BC.rd[,,i,j]))[,3]
                                       })
                           colnames(outi)=paste0("rand",1:ncol(outi))
                           data.frame(samp2.name,outi,stringsAsFactors = FALSE)
                         })
    }
    output=c(output,list(rand=BC.randm))
  }
  output
}