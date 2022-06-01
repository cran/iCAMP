qpen.cm<-function(comm,pd=NULL,pd.big.wd=NULL,pd.big.spname=NULL,tree=NULL,
                  meta.group=NULL,meta.com=NULL,meta.frequency=NULL,meta.ab=NULL,
                  ab.weight=TRUE,exclude.conspecifics=FALSE,rand.time=1000,
                  sig.bNTI=1.96,sig.rc=0.95,nworker=4,memory.G=50,project=NA,wd=getwd(),
                  output.detail=FALSE,save.bNTIRC=FALSE,
                  taxo.metric="bray", transform.method=NULL, logbase=2,
                  dirichlet=FALSE)
{
  ##########################
  # 1 # basic check and loading
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
  
  # load packages
  requireNamespace("vegan")
  requireNamespace("parallel")
  
  # creat folder if necessary
  begin=Sys.time()
  gc()
  begin.size=utils::memory.size()
  if(!is.null(pd.big.wd)){if(!dir.exists(pd.big.wd)){dir.create(path = pd.big.wd)}}
  
  if(is.na(project)){project=format(Sys.time(),format = "%Y%m%d%H%M")}
  
  ##########################
  # 2 # ID match
  
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
  if(is.null(pd))
  {
    if(is.null(tree))
    {
      stop("Neither PD nor tree was found. Calculation stopped.",date())
    }else{
      spc=iCAMP::match.name(name.check = tree$tip.label,cn.list = list(comm=comm),tree.list = list(tree=tree))
      comm=spc$comm
      tree=spc$tree
      pd.spname=tree$tip.label
    }
  }else{
    if(is.null(pd.big.wd))
    {
      colnames(pd)=rownames(pd) # sometimes R may add "X" automatically to colnames
      sp.check=match.name(cn.list = list(comm=comm),both.list = list(pd=pd))
      comm=sp.check$comm
      pd.spname=rownames(pd)
    }else{
      if(length(pd)!=1|(!file.exists(paste0(pd.big.wd,"/",pd[[1]][[1]])))){stop("the name of pd.big file must be input.")}
      if(is.null(pd.big.spname)) stop("Need species names for the pd.big.")
      spc=iCAMP::match.name(name.check = pd.big.spname,cn.list=list(comm=comm))
      comm=spc$comm
      pd.spname=pd.big.spname
    }
  }
  comms=lapply(meta.lev,function(mi){sampi=rownames(meta.group)[which(meta.group[,1]==mi)];comi=comm[which(rownames(comm) %in% sampi),,drop=FALSE];comi[,colSums(comi)>0,drop=FALSE]})
  names(comms)=meta.lev
  ##########################
  if((!is.null(meta.com)) & ((is.null(meta.frequency)) | (is.null(meta.ab))))
  {
    if(!is.list(meta.com))
    {
      if(length(meta.lev)>1){stop("meta.group has more than one metacommunities, but meta.com only defines one.")}
      meta.com=list(Meta=meta.com)
    }
    if(sum(!(meta.lev %in% names(meta.com)))>0)
    {
      stop('meta.com names must be the same as metacommunity names in meta.group.')
    }else{
      meta.com=meta.com[match(meta.lev,names(meta.com))]
    }
    if(sum(sapply(1:length(meta.lev),function(i){sum(!(colnames(comms[[i]]) %in% colnames(meta.com[[i]])))}))>0)
    {
      stop('comm has some species not included the meta.com of its metacommunity.')
    }
  }
  
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
    if(is.null(meta.com))
    {
      meta.frequency=matrix(0,nrow = length(meta.lev),ncol = ncol(comm))
      rownames(meta.frequency)=meta.lev
      colnames(meta.frequency)=colnames(comm)
      for(i in 1:length(comms))
      {
        meta.frequency[i,match(colnames(comms[[i]]),colnames(comm))]=colSums(comms[[i]]>0)
      }
    }else{
      spall=unique(unlist(lapply(meta.com,colnames)))
      meta.frequency=t(sapply(1:length(meta.com),
                            function(i)
                            {
                              outi=rep(0,length(spall))
                              outi[match(colnames(meta.com[[i]]),spall)]=colSums(meta.com[[i]]>0)
                            }))
      rownames(meta.frequency)=meta.lev
      colnames(meta.frequency)=spall
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
    if(is.null(meta.com))
    {
      meta.ab=matrix(0,nrow = length(meta.lev),ncol = ncol(comm))
      rownames(meta.ab)=meta.lev
      colnames(meta.ab)=colnames(comm)
      for(i in 1:length(comms))
      {
        meta.ab[i,match(colnames(comms[[i]]),colnames(comm))]=colMeans(comms[[i]]/rowSums(comms[[i]]))
      }
    }else{
      spall=unique(unlist(lapply(meta.com,colnames)))
      meta.ab=t(sapply(1:length(meta.com),
                              function(i)
                              {
                                outi=rep(0,length(spall))
                                outi[match(colnames(meta.com[[i]]),spall)]=colMeans(meta.com[[i]]/rowSums(meta.com[[i]]))
                              }))
      rownames(meta.ab)=meta.lev
      colnames(meta.ab)=spall
    }
  }
  meta.spool=lapply(1:nrow(meta.frequency),function(i){colnames(meta.frequency)[which(meta.frequency[i,]>0)]})
  names(meta.spool)=rownames(meta.frequency)
  
  ###############################
  # 3 # calculate phylogenetic distance
  
  if(is.null(pd))
  {
      message("---Now calculating PD ----",date())
      if(is.null(pd.big.wd))
      {
        pd.big.wd=paste0(wd,"/",project,".pdbig")
        if(!dir.exists(pd.big.wd)){dir.create(pd.big.wd)}
      }
      pdc=iCAMP::pdist.big(tree=tree,wd=pd.big.wd,nworker = nworker)
      pd=pdc$pd.file
      pd.big.spname=pdc$tip.label
    
    gc()
  }
  
  ###############################
  # 4 # calculate betaNTI
  
  if(is.null(pd.big.wd))
  {
    bNTI.M=iCAMP::bNTI.cm(comm=comm, dis=pd, nworker = nworker, memo.size.GB = memory.G,
                          meta.group=meta.group,meta.spool=meta.spool,
                          weighted = ab.weight, exclude.consp = exclude.conspecifics,
                          rand = rand.time, output.bMNTD = TRUE,sig.index="SES",
                          unit.sum = NULL, correct.special = FALSE, detail.null=TRUE)
    bNTI=bNTI.M$index
    bMNTD=bNTI.M$betaMNTD.obs
    bMNTD.rand=bNTI.M$rand
    gc()
  }else{
    bNTI.M=iCAMP::bNTI.big.cm(comm=comm, meta.group=meta.group, meta.spool = meta.spool,
                              pd.desc=pd,pd.spname=pd.big.spname,pd.wd=pd.big.wd, spname.check=TRUE,
                              nworker=nworker, memo.size.GB=memory.G, weighted=ab.weight,
                              exclude.consp=exclude.conspecifics,rand=rand.time,output.dtail=TRUE,
                              RC=FALSE, trace=TRUE)
    bNTI=bNTI.M$bNTI
    bMNTD=bNTI.M$bMNTD
    bMNTD.rand=bNTI.M$bMNTD.rand
    gc()
  }
  
  ###############################
  # 4 # calculate RC
  rc.m=iCAMP::RC.cm(comm=comm,rand=rand.time,na.zero=TRUE,nworker=nworker,
                    meta.group=meta.group,meta.frequency=meta.frequency,meta.ab=meta.ab,
                    memory.G=memory.G,weighted=ab.weight,unit.sum=NULL,
                    sig.index="RC",detail.null=TRUE,output.bray=TRUE,silent=FALSE,
                    taxo.metric=taxo.metric, transform.method=transform.method,
                    logbase=logbase,dirichlet=dirichlet)
  BC=rc.m$BC.obs
  RC=rc.m$index
  BC.rand=rc.m$rand
  
  sampc=iCAMP::match.name(both.list = list(bNTI=bNTI,RC=RC))
  bNTI=sampc$bNTI
  rcc=sampc$RC
  gc()
  
  # Proportion of different processes
  bNTI[is.na(bNTI)]=0
  bNTI.v<-as.vector(stats::as.dist(bNTI))
  rc.v<-as.vector(stats::as.dist(rcc))
  id.selectna<-(bNTI.v<=sig.bNTI&bNTI.v>=(-sig.bNTI))
  num.pair<-length(bNTI.v)
  select.h<-sum(bNTI.v>sig.bNTI)/num.pair
  select.l<-sum(bNTI.v<(-sig.bNTI))/num.pair
  disper.h<-sum(rc.v[id.selectna]>sig.rc)/num.pair
  disper.l<-sum(rc.v[id.selectna]<(-sig.rc))/num.pair
  drift<-sum(rc.v[id.selectna]<=sig.rc&rc.v[id.selectna]>=(-sig.rc))/num.pair
  res=data.frame(Heterogeneous.Selection=select.h,
                 Homogeneous.Selection=select.l,
                 Dispersal.Limitation=disper.h,
                 Homogenizing.Dispersal=disper.l,
                 Undominated=drift,num.pair)
  
  bNTI3=dist.3col(bNTI)
  RC3=dist.3col(rcc)
  bNTI.RC=data.frame(bNTI3[,3],RCv=RC3[,3]);colnames(bNTI.RC)=c("bNTI","RC")
  if(!is.null(BC)){BC3=dist.3col(BC);bNTI.RC=data.frame(BC=BC3[,3],bNTI.RC)}
  if(!is.null(bMNTD)){bMNTD3=dist.3col(bMNTD);bNTI.RC=data.frame(bMNTD=bMNTD3[,3],bNTI.RC)}
  bNTI.RC=data.frame(sample1=RC3[,1],sample2=RC3[,2],bNTI.RC,process=rep(NA,nrow(bNTI.RC)))
  id.vse=(bNTI3[,3]>sig.bNTI)
  bNTI.RC[id.vse,ncol(bNTI.RC)]="Heterogeneous.Selection"
  id.hse=(bNTI3[,3]<(-sig.bNTI))
  bNTI.RC[id.hse,ncol(bNTI.RC)]="Homogeneous.Selection"
  id.sena=(bNTI3[,3]<=sig.bNTI&bNTI3[,3]>=(-sig.bNTI))
  id.dl=(id.sena&(RC3[,3]>sig.rc))
  bNTI.RC[id.dl,ncol(bNTI.RC)]="Dispersal.Limitation"
  id.dh=(id.sena&(RC3[,3]<(-sig.rc)))
  bNTI.RC[id.dh,ncol(bNTI.RC)]="Homogenizing.Dispersal"
  id.df=(id.sena&(RC3[,3]>=(-sig.rc)&RC3[,3]<=sig.rc))
  bNTI.RC[id.df,ncol(bNTI.RC)]="Undominated"
  if(save.bNTIRC){utils::write.csv(bNTI.RC,file=paste0(wd,"/",project,".QPEN.detail.csv"))}
  end=Sys.time()
  end.size=utils::memory.size()
  setting=data.frame(project=project,ab.weighted=ab.weight,
                     rand.time=rand.time,exclude.conspecifics=exclude.conspecifics,
                     nworker=nworker,memory.set=paste(memory.G,"Gb",sep=" "),
                     memory.use=paste((end.size-begin.size),"Mb",sep=" "),time.use=(end-begin))
  if(!output.detail){output=list(ratio=res,result=bNTI.RC)}else{output=list(ratio=res,result=bNTI.RC,pd=pd,bMNTD=bMNTD,BC=BC,bMNTD.rand=bMNTD.rand,BC.rand=BC.rand,setting=setting)}
  output
}