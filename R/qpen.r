qpen<-function(comm=NULL,pd=NULL,pd.big.wd=NULL,pd.big.spname=NULL,tree=NULL,bNTI=NULL,RC=NULL,
               ab.weight=TRUE,meta.ab=NULL,exclude.conspecifics=FALSE,rand.time=1000,
               sig.bNTI=1.96,sig.rc=0.95,nworker=4,memory.G=50,project=NA,wd=getwd(),
               output.detail=FALSE,save.bNTIRC=FALSE)
{
  # quantify processes according to Stegen et al 2015
  
  if(is.na(project)){project=format(Sys.time(),format = "%Y%m%d%H%M")}
  if(.Platform$OS.type=="windows")
  {
    if(utils::memory.limit()<memory.G*1024)
    {
      memotry=try(utils::memory.limit(size=memory.G*1024),silent = TRUE)
      if(class(memotry)=="try-error"){warning(memotry[1])}
    }
  }
  
  ## load packages
  requireNamespace("vegan")
  requireNamespace("parallel")
  
  ## creat folder if necessary
  begin=Sys.time()
  gc()
  begin.size=utils::memory.size()
  if(!is.null(pd.big.wd)){if(!dir.exists(pd.big.wd)){dir.create(path = pd.big.wd)}}
  if(is.null(bNTI)|is.null(RC)){if(is.null(comm)) stop("comm must be input, if either bNTI or RC is mising.")}
  if(is.null(bNTI))
  {
    if(is.null(pd))
    {
      if(is.null(tree)){stop("Neither bNTI, PD nor tree was found. Calculation stopped.",date())}else{
        spc=iCAMP::match.name(name.check = tree$tip.label,cn.list = list(comm=comm),tree.list = list(tree=tree))
        comm=spc$comm
        tree=spc$tree
        message("---Now calculating PD ----",date())
        if(is.null(pd.big.wd))
        {
          pd.big.wd=paste0(wd,"/",project,".pdbig")
          if(!dir.exists(pd.big.wd)){dir.create(pd.big.wd)}
        }
        pdc=iCAMP::pdist.big(tree=tree,wd=pd.big.wd,nworker = nworker)
        pd=pdc$pd.file
        pd.big.spname=pdc$tip.label
      }
      gc()
    }else{
      if(is.null(pd.big.wd))
      {
        colnames(pd)=rownames(pd) # sometimes R may add "X" automatically to colnames
        sp.check=match.name(cn.list = list(comm=comm),both.list = list(pd=pd))
        comm<-sp.check$comm
        pd<-sp.check$pd
      }else{
        if(length(pd)!=1|(!file.exists(paste0(pd.big.wd,"/",pd[[1]][[1]])))){stop("the name of pd.big file must be input.")}
        if(is.null(pd.big.spname)) stop("Need species names for the pd.big.")
        spc=iCAMP::match.name(name.check = pd.big.spname,cn.list=list(comm=comm))
        comm=spc$comm
      }
    }
    if(is.null(pd.big.wd))
    {
      bNTI.M=iCAMP::bNTIn.p(comm=comm, dis=pd, nworker = nworker, memo.size.GB = memory.G,
                          weighted = ab.weight, exclude.consp = exclude.conspecifics,
                          rand = rand.time, output.bMNTD = TRUE,sig.index="SES",
                          unit.sum = NULL, correct.special = FALSE, detail.null=TRUE)
      bNTI=bNTI.M$index
      bMNTD=bNTI.M$betaMNTD.obs
      bMNTD.rand=bNTI.M$rand
      gc()
    }else{
      bNTI.M=iCAMP::bNTI.big(comm=comm, meta.group=NULL, pd.desc=pd,
                             pd.spname=pd.big.spname,pd.wd=pd.big.wd, spname.check=TRUE,
                             nworker=nworker, memo.size.GB=memory.G, weighted=ab.weight,
                             exclude.consp=exclude.conspecifics,rand=rand.time,output.dtail=TRUE,
                             RC=FALSE, trace=TRUE)
      bNTI=bNTI.M$bNTI
      bMNTD=bNTI.M$bMNTD
      bMNTD.rand=bNTI.M$bMNTD.rand
      gc()
    }
  }else{
    bMNTD=NULL;bMNTD.rand=NULL
  }
  
  if(is.null(RC))
  {
    rc.m<-RC.pc(comm=comm,rand=rand.time,na.zero=TRUE,nworker=nworker,
                memory.G=memory.G,weighted=ab.weight,unit.sum=NULL,
                meta.ab=meta.ab,sig.index="RC",
                detail.null=TRUE,output.bray=TRUE)
    BC=rc.m$BC.obs
    RC=rc.m$index
    BC.rand=rc.m$rand
  }else{
    BC=NULL;BC.rand=NULL
  }
  
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