pdist.big<-function(tree,wd=getwd(),tree.asbig=FALSE,output=FALSE,
                    nworker=4,nworker.pd=nworker,memory.G=50,time.count=FALSE,
                    treepath.file="path.rda", pd.spname.file="pd.taxon.name.csv",
                    pd.backingfile="pd.bin", pd.desc.file="pd.desc",
                    tree.backingfile="treeinfo.bin", tree.desc.file="treeinfo.desc")
{
  if(time.count) t1=Sys.time()
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memory.G*1024)
    {
      memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
      if(inherits(memotry,"try-error")){warning(memotry[1])}
    }
  }
  
  requireNamespace("bigmemory")
  requireNamespace("parallel")
  
  if(file.exists(wd))
  {
    filenames=c(tree.backingfile,tree.desc.file,pd.backingfile,pd.desc.file,pd.spname.file)
    files.wd=list.files(path=wd,full.names = FALSE)
    file.error=filenames[filenames %in% files.wd]
    if(length(file.error)>0) stop("The work directory (wd) you assigned has already had some files: ",paste(file.error,collapse = ", "),". \n Please change a wd. You may need to delete the old files to save space. ",date())
  }else{dir.create(path = wd)}
  
  if(tree.asbig)
  {
    message("Now saving tree informaiton as big matrix file. ",date())
    tree.info=bigmemory::big.matrix(nrow=length(tree$edge.length),ncol=3,type = "double",
                                    backingfile = tree.backingfile,backingpath = wd,
                                    descriptorfile = tree.desc.file,shared = TRUE)
    tree.info[]=cbind(tree$edge,tree$edge.length)
    gc()
    tip.num=length(tree$tip)
    pathf1<-function(i,tr)
    {
      requireNamespace("bigmemory")
      edge<-bigmemory::attach.big.matrix(dget(tr))
      pathx<-list()
      j=1
      pathx[[1]]=edge[mwhich(edge,2,i,'eq'),1]
      pathx[[2]]=edge[mwhich(edge,2,i,'eq'),3]
      while(length(edge[mwhich(edge,2,pathx[[1]][j],'eq'),])!=0){
        pathx[[1]][j+1]=edge[mwhich(edge,2,pathx[[1]][j],'eq'),1]
        pathx[[2]][j+1]=edge[mwhich(edge,2,pathx[[1]][j],'eq'),3]
        j=j+1
      }
      pathx
    }
    
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
    message("Now computing path. begin at ", date(),". Please wait...")
    path<-parallel::parLapply(c1,1:tip.num,pathf1,tr=paste0(wd,"/",tree.desc.file))
    parallel::stopCluster(c1)
    gc()
  }else{
    edge<-tree$edge
    edge.len<-tree$edge.length
    tip.num=length(tree$tip)
    
    pathf<-function(i,edge,edge.len)
    {
      pathx<-list()
      j=1
      pathx[[1]]=edge[edge[,2]==i,1]
      pathx[[2]]=edge.len[edge[,2]==i]
      while(length(edge[edge[,2]==pathx[[1]][j],])!=0){
        pathx[[1]][j+1]=edge[edge[,2]==pathx[[1]][j],1]
        pathx[[2]][j+1]=edge.len[edge[,2]==pathx[[1]][j]]
        j=j+1
      }
      pathx
    }
    
    
    if(tip.num>(5*nworker))
    {
      t1=Sys.time()
      c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
      if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
      if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
      t2=Sys.time();tmkc=as.numeric(difftime(t2,t1,units = "hours"))
      message("Setting parallel cluster for path computing cost ",format(t2-t1),".  ",date())
      path=list()
      path[1:(5*nworker)]<-parallel::parLapply(c1,1:(5*nworker),pathf,edge=edge,edge.len=edge.len)
      t3=Sys.time();tpara=as.numeric(difftime(t3,t2,units = "hours"))
      message("Parallel for ", 5*nworker," tips cost ",format(t3-t2),". ",date())
      parallel::stopCluster(c1)
      gc()
      t4=Sys.time();tclose=as.numeric(difftime(t4,t3,units = "hours"))
      t.total=tmkc+((tpara/5)*(ceiling(tip.num/nworker)-5))+tclose
      message("Path computing by parallel may take ", t.total," hours. ",date())
      t5=Sys.time()
      c2<-try(parallel::makeCluster(nworker,type="PSOCK"))
      if(inherits(c2,"try-error")){c2 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
      if(inherits(c2,"try-error")){c2 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
      message("Now computing path for the rest ",tip.num-(5*nworker)," tips. begin at ", date(),". Please wait...")
      path[(5*nworker+1):tip.num]<-parallel::parLapply(c2,(5*nworker+1):tip.num,pathf,edge=edge,edge.len=edge.len)
      parallel::stopCluster(c2)
      gc()
      message("Computing path for the rest ",tip.num-(5*nworker)," tips actually took ", format(Sys.time()-t5),". ",date())
    }else{
      t5=Sys.time()
      c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
      if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
      if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
      message("Now computing path. begin at ", date(),". Please wait...")
      path<-parallel::parLapply(c1,1:tip.num,pathf,edge=edge,edge.len=edge.len)
      parallel::stopCluster(c1)
      gc()
      message("Computing path took ", format(Sys.time()-t5),". ",date())
    }
  }
  save(path,file = paste0(wd,"/",treepath.file))
  if(!file.exists(paste0(wd,"/",pd.desc.file)))
  {
    message("Now setting big matrix file on the disk. ",date())
    pd.big=bigmemory::big.matrix(nrow = tip.num,ncol = tip.num, type="double",
                                 backingfile = pd.backingfile, backingpath = wd,
                                 descriptorfile = pd.desc.file,shared = TRUE)
    pd.big[]=0
  }
  utils::write.csv(data.frame(taxon.name=tree$tip.label),file=paste0(wd,"/",pd.spname.file))
  
  pdf<-function(m, tip.num, path, pdsc)
  {
    message("now m=",m,". ",date())
    pdist<-bigmemory::attach.big.matrix(dget(pdsc))
    pdff<-function(n)
    {
      a=path[[m]][[1]]
      b=path[[n]][[1]]
      al=path[[m]][[2]]
      bl=path[[n]][[2]]
      sum(al[c(TRUE,is.na(match(a,b)[1:(length(a)-1)]))])+sum(bl[c(TRUE,is.na(match(b,a)[1:(length(b)-1)]))])
    }
    pdist[(m+1):tip.num,m]<-pdist[m,(m+1):tip.num]<-sapply((m+1):tip.num, pdff)
    invisible()
  }
  
  if((tip.num-1)>(5*nworker.pd))
  {
    t1=Sys.time()
    c3<-try(parallel::makeCluster(nworker.pd,type="PSOCK"))
    if(inherits(c3,"try-error")){c3 <- try(parallel::makeCluster(nworker.pd, setup_timeout = 0.5))}
    if(inherits(c3,"try-error")){c3 <- parallel::makeCluster(nworker.pd, setup_strategy = "sequential")}
    t2=Sys.time();tmkc=as.numeric(difftime(t2,t1,units = "hours"))
    message("Setting parallel cluster for pdist computing cost ",format(t2-t1),".  ",date())
    pd.cal=list()
    pd.cal[1:(5*nworker.pd)]<-parallel::parLapply(c3,1:(5*nworker.pd),pdf,tip.num=tip.num,path=path, pdsc=paste0(wd,"/",pd.desc.file))
    t3=Sys.time();tpara=as.numeric(difftime(t3,t2,units = "hours"))
    message("Parallel computing Pdist for the first ",5*nworker.pd ," runs cost ",format(t3-t2),". ",date())
    parallel::stopCluster(c3)
    gc()
    t4=Sys.time();tclose=as.numeric(difftime(t4,t3,units = "hours"))
    tfactor=(1-((((ceiling((tip.num-1)/nworker.pd))-1)*nworker.pd+1)/(2*(tip.num-1))))*ceiling((tip.num-1)/nworker.pd)-5
    t.total=tmkc+((tpara/5)*tfactor)+tclose
    message("The rest Pdist computing by parallel may take ", t.total," hours. ",date())
    t5=Sys.time()
    c4<-try(parallel::makeCluster(nworker.pd,type="PSOCK"))
    if(inherits(c4,"try-error")){c4 <- try(parallel::makeCluster(nworker.pd, setup_timeout = 0.5))}
    if(inherits(c4,"try-error")){c4 <- parallel::makeCluster(nworker.pd, setup_strategy = "sequential")}
    pd.cal[(5*nworker.pd+1):(tip.num-1)]<-parallel::parLapply(c4,(5*nworker.pd+1):(tip.num-1),pdf,tip.num=tip.num,path=path, pdsc=paste0(wd,"/",pd.desc.file))
    parallel::stopCluster(c4)
    gc()
    message("Computing pdist for the rest ",tip.num-(5*nworker.pd)," tips actually took ", format(Sys.time()-t5),". ",date())
  }else{
    t5=Sys.time()
    c3<-try(parallel::makeCluster(nworker.pd,type="PSOCK"))
    if(inherits(c3,"try-error")){c3 <- try(parallel::makeCluster(nworker.pd, setup_timeout = 0.5))}
    if(inherits(c3,"try-error")){c3 <- parallel::makeCluster(nworker.pd, setup_strategy = "sequential")}
    message("Now computing pdist. begin at ", date(),". Please wait...")
    pd.cal<-parallel::parLapply(c3,1:(tip.num-1),pdf,tip.num=tip.num,path=path, pdsc=paste0(wd,"/",pd.desc.file))
    parallel::stopCluster(c3)
    message("Computing pdist actually took ", format(Sys.time()-t5),". ",date())
    gc()
  }
  
  if(output)
  {
    pdist<-bigmemory::attach.big.matrix(dget(paste0(wd,"/",pd.desc.file)))
    output=pdist[]
    rownames(output)=tree$tip.label
    colnames(output)=tree$tip.label
  }else{
    output=list(tip.label=tree$tip.label,pd.wd=wd,pd.file=pd.desc.file,pd.name.file=pd.spname.file)
  }
  
  if(time.count){t2=format(Sys.time()-t1);message("----Phylogenetic distance calculation costed ",t2," in total----")}
  output
}