dniche<-function(env,comm,method=c("ab.overlap","niche.value","prefer.overlap"),
                 nworker=4,memory.G=50,out.dist=FALSE,bigmemo=TRUE,nd.wd=getwd(),
                 nd.spname.file="nd.names.csv",detail.file="ND.res.rda")
{
  #source(paste(code.wd,"/match.name.r",sep = ""))
  checksamp=match.name(name.check = rownames(env),rn.list = list(env=env,comm=comm))
  env=checksamp$env
  comm=checksamp$comm
  envnum=ncol(env)
  spname=colnames(comm)
  if(is.null(colnames(env))){colnames(env)=paste0("env",1:ncol(env))}
  if(is.null(nd.wd)){nd.wd=getwd()}
  
  res=list()
  res$bigmemo=bigmemo
  res$nd=list()
  if(bigmemo)
  {
    requireNamespace("bigmemory")
    res$nd.wd=nd.wd
    res$names=spname
    utils::write.csv(data.frame(names=spname),file = paste0(nd.wd,"/",nd.spname.file))
    ndbig=list()
    for(j in 1:ncol(env))
    {
      ndbig[[j]] = bigmemory::big.matrix(nrow = length(spname), ncol = length(spname), type = "double",
                                         backingfile = paste0(colnames(env)[j],".ND.bin"),
                                         backingpath = nd.wd,
                                         descriptorfile = paste0(colnames(env)[j],".ND.desc"),
                                         shared = TRUE)
      ndbig[[j]][]=0
      res$nd[[j]]=paste0(colnames(env)[j],".ND.desc")
    }
  }
  
  if(method[1]=="niche.value")
  {
    comts<-(t(comm)/colSums(comm))
    nv<-(as.matrix(comts) %*% as.matrix(env))
    for(j in 1:ncol(env))
    {
      ndjm=as.matrix(stats::dist(nv[,j]))
      if(bigmemo)
      {
        ndbig[[j]][] = ndjm/max(ndjm)
      }else{
        res$nd[[j]]=ndjm/max(ndjm)
      }
    }
  }else if(method[1] %in% c("ab.overlap","prefer.overlap")){
    if(method[1]=="ab.overlap")
    {
      comp=comm
    }else{
      comp=t(t(comm)/colSums(comm))
    }
    dens<-function(i,envj,comp)
    {
      stats::density(envj,weights = comp[,i],from=min(envj),to=max(envj))$y
    }
    dio<-function(i,den,res,j)
    {
      den.b<-den.max<-den[,(i+1):ncol(den),drop=FALSE]
      den.a<-den.am<-matrix(den[,i],nrow=nrow(den),ncol=ncol(den.b))
      den.max[den.max<den.a]=0
      den.am[den.am<=den.b]=0
      den.max=den.max+den.am
      den.dif=abs(den.b-den.a)
      ndio=colSums(den.dif)/colSums(den.max)
      if(res$bigmemo)
      {
        requireNamespace("bigmemory")
        ndj=bigmemory::attach.big.matrix(dget(res$nd[[j]]))
        ndj[i,i]=0
        ndj[i,(i+1):ncol(den)]=ndio
        ndj[(i+1):ncol(den),i]=ndio
        out=i
      }else{
        out=rep(0,ncol(den))
        out[(i+1):ncol(den)]=ndio
      }
      out
    }
    
    requireNamespace("parallel")
    if(.Platform$OS.type=="windows")
    {
      if(utils::memory.limit()<memory.G*1024)
      {
        memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
        if(class(memotry)=="try-error"){warning(memotry[1])}
      }
    }
    
    for(j in 1:envnum)
    {
      dens1=stats::density(env[,j],weights = comp[,1],from=min(env[,j]),to=max(env[,j]))
      den1=data.frame(dens1$y)
      rownames(den1)<-dens1$x
      if(nworker==1)
      {
        message("Now computing density model. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        den<-lapply(2:ncol(comm),dens,envj=env[,j],comp=comp)
      }else{
        c1<-parallel::makeCluster(nworker,type="PSOCK")
        message("Now parallel computing density model. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        den<-parallel::parLapply(c1,2:ncol(comm),dens,envj=env[,j],comp=comp)
        parallel::stopCluster(c1)
      }
      
      den=matrix(unlist(den),nrow=length(den[[1]]))
      den=as.matrix(cbind(den1,den))
      gc()
      
      if(nworker==1)
      {
        message("Now computing niche distance. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        dis<-lapply(c0,1:(ncol(den)-1),dio,den=den,res=res,j=j)
      }else{
        c0<-parallel::makeCluster(nworker,type="PSOCK")
        message("Now parallel computing niche distance. j=",j," in ",envnum,". begin at ", date(),". Please wait...")
        dis<-parallel::parLapply(c0,1:(ncol(den)-1),dio,den=den,res=res,j=j)
        parallel::stopCluster(c0)
      }
      if(!bigmemo)
      {
        dis=data.frame(matrix(unlist(dis),nrow=length(dis[[1]])))
        dis=cbind(dis,rep(0,ncol(comm)))
        dis=dis+t(dis)
        colnames(dis)<-rownames(dis)<-spname
        if(out.dist){dis=stats::as.dist(dis)}
        res$nd[[j]]=dis
      }
    }
  }
  names(res$nd)=colnames(env)
  res$method=method[1]
  if(bigmemo){if(!is.null(detail.file)){save(res,file = paste0(nd.wd,"/",detail.file))}}
  res
}