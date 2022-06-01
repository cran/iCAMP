icamp.big<-function(comm,tree,pd.desc=NULL, pd.spname=NULL, pd.wd=getwd(),
                    rand=1000,prefix="iCAMP",ds=0.2,pd.cut=NA,sp.check=TRUE,
                    phylo.rand.scale=c("within.bin","across.all","both"),
                    taxa.rand.scale=c("across.all","within.bin","both"),
                    phylo.metric=c("bMPD","bMNTD","both", "bNRI", "bNTI"),
                    sig.index=c("Confidence","SES.RC","SES","RC"),
                    bin.size.limit=24,nworker=4,memory.G=50,
                    rtree.save=FALSE,detail.save=TRUE,qp.save=TRUE,detail.null=FALSE,
                    ignore.zero=TRUE,output.wd=getwd(),
                    correct.special=TRUE,unit.sum=rowSums(comm),
                    special.method=c("depend","MPD","MNTD","both"),
                    ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975,
                    omit.option=c("no","test","omit"),meta.ab=NULL,
                    treepath.file="path.rda", pd.spname.file="pd.taxon.name.csv",
                    pd.backingfile="pd.bin", pd.desc.file="pd.desc",
                    taxo.metric="bray", transform.method=NULL,
                    logbase=2, dirichlet=FALSE, d.cut.method=c("maxpd","maxdroot"))
{
  # v20200725 add sig.index, conf.cut, detail.null
  # v20210417 add taxo.metric, transform.method, logbase, dirichlet
  comm.old=comm
  comm=as.matrix(comm.old)
  if(sum(comm!=comm.old)>0){stop("Strangely, comm changed after using as.matrix. Please check the comm.")}
  if(max(rowSums(comm,na.rm = TRUE))<=1 & (!dirichlet))
  {
    warning("The values in comm are less than 1, thus considered as proportional data, Dirichlet distribution is used to assign abundance in null model.")
    dirichlet=TRUE
  }
  taxo.name=paste0(toupper(substr(taxo.metric,1,1)),substring(taxo.metric,2))
  
  t1=Sys.time()
  phylo.rand.scale=phylo.rand.scale[1]
  taxa.rand.scale=taxa.rand.scale[1]
  phylo.metric=phylo.metric[1]
  sig.index=sig.index[1]
  special.method=special.method[1]
  omit.option=omit.option[1]
  if(!(sig.index %in% c("SES.RC","Confidence","RC","SES"))){stop("wrong sig.index for icamp.big.")}
  
  requireNamespace("vegan")
  requireNamespace("ape")
  requireNamespace("bigmemory")
  if(detail.null){detail.rand=list();detail.special=list()}
  
  if(!dir.exists(pd.wd)){dir.create(pd.wd)}
  if(is.null(pd.desc))
  {
    if(file.exists(paste0(pd.wd,"/",pd.desc.file)))
    {
      warning("Attention: the pd.wd already has a pd.desc.file, which is directly used. Please double check whether this is the phylogenetic distance file you need!")
      pd.big=list()
      pd.big$tip.label=utils::read.csv(paste0(pd.wd,"/",pd.spname.file),row.names = 1,stringsAsFactors = FALSE)[,1]
      pd.big$pd.wd=pd.wd
      pd.big$pd.file=pd.desc.file
      pd.big$pd.name.file=pd.spname.file
    }else{
      pd.big=iCAMP::pdist.big(tree = tree, wd=pd.wd, nworker = nworker,
                              treepath.file=treepath.file, pd.spname.file=pd.spname.file,
                              pd.backingfile=pd.backingfile, pd.desc.file=pd.desc.file)
    }
    pd.desc=pd.big$pd.file
    pd.spname=pd.big$tip.label
    pd.wd=pd.big$pd.wd
  }
  
  pd=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
  diagid=matrix(c(1:nrow(pd),1:nrow(pd)),nrow=nrow(pd),ncol=2)
  pd[diagid]=0
  
  if(!is.rooted(tree))
  {
    tree.rt=iCAMP::midpoint.root.big(tree, pd.desc=pd.desc, pd.spname=pd.spname,
                                     pd.wd=pd.wd, nworker=nworker)
    tree=tree.rt$tree
    maxpd=tree.rt$max.pd
    if(rtree.save) ape::write.tree(tree,file = paste0(output.wd,"/",prefix,".tree.rooted.nwk"))
  }else{maxpd=NA}
  
  if(ignore.zero) comm=comm[rowSums(comm)>0,colSums(comm)>0]
  
  if(sp.check)
  {
    checksp=iCAMP::match.name(name.check = pd.spname,cn.list = list(comm=comm),tree.list=list(tree=tree))
    comm=checksp$comm
    tree=checksp$tree
  }
  samp.num=nrow(comm)
  
  if(is.na(pd.cut)){d.cut=NULL}else{d.cut=pd.cut}
  
  message("----------Now binning-----------------",date())
  taxabin3=iCAMP::taxa.binphy.big(tree,pd.desc=pd.desc, pd.spname=pd.spname,
                                  pd.wd=pd.wd,outgroup.tip=NA,outgroup.rm=TRUE,
                                  d.cut=d.cut,ds=ds,bin.size.limit=bin.size.limit,
                                  nworker=nworker,d.cut.method=d.cut.method)
  gc()
  if(detail.save)
  {
    detail=list()
    detail$taxabin=taxabin3
  }
  if(omit.option[1]=="no")
  {
    message("----------Now binning com and sp without omitting small bins-----------------",date())
    sp.bin=taxabin3$sp.bin[,3,drop=FALSE]
  }else if(omit.option[1] %in% c("test","omit")){
    message("----------Now test what will be lost if omitting small bins-----------------",date())
    str.bin=taxabin3$state.strict
    rm.bins=str.bin$bin.strict.id[which(str.bin$bin.strict.taxa.num<bin.size.limit)]
    if(length(rm.bins)==0){message("No bin is smaller than limitation. Good to go!")}else{
      rm.sps=taxabin3$sp.bin[which(taxabin3$sp.bin$bin.id.strict %in% rm.bins),1,drop=FALSE]
      com.rm=comm[,match(rownames(rm.sps),colnames(comm)),drop=FALSE]
      rm.binnum=paste0(length(rm.bins), " in ",length(taxabin3$bin.strict.sp)," (",round(length(rm.bins)/length(taxabin3$bin.strict.sp),4)*100,"%)")
      com.rmra=com.rm/rowSums(comm)
      rm.absum=sum(colMeans(com.rmra))
      sum.line=c(rm.binnum,rm.absum,rowSums(com.rmra))
      omit.out=rbind(Summary=sum.line,
                     data.frame(rm.sps,Relative.abundance=colMeans(com.rmra),t(com.rmra),stringsAsFactors = FALSE))
      message("Omit bin number ",rm.binnum,", with total relative abundance ",round(rm.absum,4)*100,"%.")
      if(detail.save)
      {
        detail$omit.bins=omit.out
      }
      utils::write.csv(omit.out,file=paste0(output.wd,"/",prefix,".OmitBinInform.csv"))
    }
    if(omit.option[1]=="test")
    {
      message("Since omit.option is test, iCAMP is stopped. please check OmitBinInform.csv in output folder.")
      return(omit.out)
    }else{
      if(rm.absum>0.95){stop("The omitted bins have too high relative abundance, ",round(rm.absum,4)*100,"%, iCAMP is stopped.")}
      if(rm.absum>0.2){warning("The omitted bins have a total relative abundance of ",round(rm.absum,4)*100,"%.")}
      message("----------Small bins are omitted-----------------",date())
      sp.bin=taxabin3$sp.bin[which(!(taxabin3$sp.bin$bin.id.strict %in% rm.bins)),1,drop=FALSE]
      comm=comm[,match(rownames(sp.bin),colnames(comm)),drop=FALSE]
    }
  }else{stop("omit.option input error.")}
  
  bin.lev=levels(as.factor(sp.bin[,1]))
  bin.num=length(bin.lev)
  com.bin=lapply(1:bin.num, function(i){comm[,match(rownames(sp.bin)[which(sp.bin==bin.lev[i])],colnames(comm))]})
  # thanks to adityabandla, I revised sp.bin==i to sp.bin==bin.lev[i], otherwise omit.option='omit' will lead to wrong results.
  # see https://github.com/DaliangNing/iCAMP1/issues/9 for detailed discussion.
  gc()
  pdid.bin=lapply(1:bin.num, function(i){match(rownames(sp.bin)[which(sp.bin==bin.lev[i])],pd.spname)})
  gc()
  bin.abund=sapply(1:bin.num,function(i){rowSums(com.bin[[i]])})
  colnames(bin.abund)=paste("bin",bin.lev,sep = "")
  
  if(special.method=="depend")
  {
    special.method.bNRI="MPD"
    special.method.bNTI="MNTD"
  }else{
    special.method.bNRI<-special.method.bNTI<-special.method
  }
  
  if(sig.index=="SES.RC")
  {
    sig.ind.phy="SES";sig.ind.tax="RC"
  }else if(sig.index=="SES"){
    sig.ind.phy<-sig.ind.tax<-"SES"
  }else if(sig.index=="RC"){
    sig.ind.phy<-sig.ind.tax<-"RC"
  }else if(sig.index=="Confidence"){
    sig.ind.phy<-sig.ind.tax<-"Confidence"
  }else{stop("wrong sig.index in icamp.big.")}
  
  if(phylo.rand.scale %in% c("within.bin","both"))
  {
    if(phylo.metric %in% c("bMPD","bNRI","both"))
    {
      bNRI1res=lapply(1:bin.num,
                      function(i)
                      {
                        message("----bMPD bin i=",i," in ",bin.num," ---- ",date())
                        pd.bini=pd[pdid.bin[[i]],pdid.bin[[i]]]
                        rownames(pd.bini)<-colnames(pd.bini)<-colnames(com.bin[[i]])
                        iCAMP::bNRIn.p(comm = com.bin[[i]], dis = pd.bini,nworker=nworker,
                                       memo.size.GB = memory.G,weighted = TRUE,rand = rand,
                                       sig.index=sig.ind.phy,unit.sum = unit.sum,output.bMPD = detail.null,
                                       correct.special = correct.special,detail.null=detail.null,
                                       special.method=special.method.bNRI,ses.cut=ses.cut,
                                       rc.cut=rc.cut,conf.cut=conf.cut,dirichlet=dirichlet)
                      })
      gc()
      bNRI1=lapply(1:length(bNRI1res),function(i){bNRI1res[[i]]$index})
      if(sig.ind.phy=="Confidence"){objname="CbMPDi"}else if(sig.ind.phy=="SES"){objname="bNRIi"}else if(sig.ind.phy=="RC"){objname="RCbMPDi"}
      bNRI1.m=iCAMP::dist.bin.3col(bNRI1,obj.name = objname)
      if(detail.save){detail$SigbMPDi=bNRI1.m}
      if(detail.null)
      {
        detail.rand$bMPD.obs=lapply(1:length(bNRI1res),function(i){bNRI1res[[i]]$betaMPD.obs})
        detail.rand$bMPDi.rand=lapply(1:length(bNRI1res),function(i){bNRI1res[[i]]$rand})
        
        if(correct.special)
        {
          special.bNRI1=lapply(1:length(bNRI1res[[1]]$special.crct),
                               function(i)
                               {
                                 iCAMP::dist.bin.3col(lapply(1:length(bNRI1res),function(j){bNRI1res[[j]]$special.crct[[i]]}))
                               })
          names(special.bNRI1)=names(bNRI1res[[1]]$special.crct)
        }else{special.bNRI1=NULL}
        detail.special$SigbMPDi=special.bNRI1
      }
    }
    if(phylo.metric %in% c("bMNTD","bNTI","both"))
    {
      bNTI1res=lapply(1:bin.num,
                      function(i)
                      {
                        message("----bMNTD bin i=",i," in ",bin.num," ---- ",date())
                        pd.bini=pd[pdid.bin[[i]],pdid.bin[[i]]]
                        rownames(pd.bini)<-colnames(pd.bini)<-colnames(com.bin[[i]])
                        iCAMP::bNTIn.p(comm=com.bin[[i]],dis=pd.bini,nworker = nworker,
                                       memo.size.GB = memory.G,weighted = TRUE,rand = rand,
                                       sig.index=sig.ind.phy,unit.sum=unit.sum,output.bMNTD = detail.null,
                                       correct.special = correct.special,detail.null=detail.null,
                                       special.method=special.method.bNTI,ses.cut=ses.cut,
                                       rc.cut=rc.cut,conf.cut=conf.cut,dirichlet=dirichlet)
                      })
      gc()
      bNTI1=lapply(1:length(bNTI1res),function(i){bNTI1res[[i]]$index})
      if(sig.ind.phy=="Confidence"){objname="CbMNTDi"}else if(sig.ind.phy=="SES"){objname="bNTIi"}else if(sig.ind.phy=="RC"){objname="RCbMNTDi"}
      bNTI1.m=iCAMP::dist.bin.3col(bNTI1, obj.name = objname)
      if(detail.save){detail$SigbMNTDi=bNTI1.m}
      if(detail.null)
      {
        detail.rand$bMNTD.obs=lapply(1:length(bNTI1res),function(i){bNTI1res[[i]]$betaMNTD.obs})
        detail.rand$bMNTDi.rand=lapply(1:length(bNTI1res),function(i){bNTI1res[[i]]$rand})
        if(correct.special)
        {
          special.bNTI1=lapply(1:length(bNTI1res[[1]]$special.crct),
                               function(i)
                               {
                                 iCAMP::dist.bin.3col(lapply(1:length(bNTI1res),function(j){bNTI1res[[j]]$special.crct[[i]]}))
                               })
          names(special.bNTI1)=names(bNTI1res[[1]]$special.crct)
        }else{special.bNTI1=NULL}
        detail.special$SigbMNTDi=special.bNTI1
      }
    }
  }
  if(phylo.rand.scale %in% c("across.all","both"))
  {
    if(phylo.metric %in% c("bMPD","bNRI","both"))
    {
      bNRI2res=iCAMP::bNRI.bin.big(comm = comm,pd.desc=pd.desc, pd.spname=pd.spname, pd.wd=pd.wd,
                                   pdid.bin=pdid.bin, sp.bin = sp.bin,spname.check = FALSE,
                                   nworker = nworker,memo.size.GB = memory.G,weighted = TRUE,
                                   rand = rand,output.bMPD=detail.null,sig.index=sig.ind.phy,
                                   unit.sum=unit.sum,correct.special = correct.special,
                                   detail.null=detail.null,special.method=special.method.bNRI,
                                   ses.cut=ses.cut,rc.cut=rc.cut,conf.cut=conf.cut,dirichlet=dirichlet)
      if(sig.ind.phy=="Confidence"){objname="CbMPDa"}else if(sig.ind.phy=="SES"){objname="bNRIa"}else if(sig.ind.phy=="RC"){objname="RCbMPDa"}
      bNRI2.m=iCAMP::dist.bin.3col(bNRI2res$index, obj.name = objname)
      if(detail.save){detail$SigbMPDa=bNRI2.m}
      if(detail.null)
      {
        detail.rand$bMPD.obs=bNRI2res$betaMPD.obs
        detail.rand$bMPDa.rand=bNRI2res$rand
        detail.special$SigbMPDa=bNRI2res$special.crct
      }
    }
    
    if(phylo.metric %in% c("bMNTD","bNTI","both"))
    {
      bNTI2res=iCAMP::bNTI.bin.big(comm = comm, pd.desc=pd.desc, pd.spname=pd.spname, pd.wd=pd.wd,
                                   pdid.bin=pdid.bin, sp.bin = sp.bin,spname.check = FALSE,
                                   nworker = nworker,memo.size.GB = memory.G,weighted = TRUE,
                                   rand = rand,output.bMNTD = detail.null,sig.index=sig.ind.phy,
                                   unit.sum=unit.sum,correct.special = correct.special,
                                   detail.null=detail.null,special.method=special.method.bNTI,
                                   ses.cut=ses.cut,rc.cut=rc.cut,conf.cut=conf.cut,dirichlet=dirichlet)
      if(sig.ind.phy=="Confidence"){objname="CbMNTDa"}else if(sig.ind.phy=="SES"){objname="bNTIa"}else if(sig.ind.phy=="RC"){objname="RCbMNTDa"}
      bNTI2.m=iCAMP::dist.bin.3col(bNTI2res$index, obj.name = objname)
      if(detail.save){detail$SigbMNTDa=bNTI2.m}
      if(detail.null)
      {
        detail.rand$bMNTD.obs=bNTI2res$betaMNTD.obs
        detail.rand$bMNTDa.rand=bNTI2res$rand
        detail.special$SigbMNTDa=bNTI2res$special.crct
      }
    }
  }
  
  ###################################################################
  #####################################################################  
  if(taxa.rand.scale %in% c("within.bin","both"))
  {
    RC1res=lapply(1:bin.num,
                  function(i)
                  {
                    message("----Bray bin i=",i," in ",bin.num," ---- ",date())
                    iCAMP::RC.pc(comm = com.bin[[i]],rand = rand,
                                 nworker = nworker,memory.G = memory.G,
                                 unit.sum = unit.sum,meta.ab = meta.ab,
                                 sig.index=sig.ind.tax,detail.null=detail.null,
                                 output.bray=detail.null,taxo.metric=taxo.metric,
                                 transform.method=transform.method, logbase=logbase,
                                 dirichlet=dirichlet)
                  })
    RC1=lapply(1:length(RC1res),function(i){RC1res[[i]]$index})
    if(sig.ind.tax=="Confidence"){objname=paste0("C",taxo.metric,"i")}else if(sig.ind.tax=="SES"){objname=paste0("SES",taxo.metric,"i")}else if(sig.ind.tax=="RC"){objname=paste0("RC",taxo.metric,"i")}
    RC1.m=iCAMP::dist.bin.3col(RC1, obj.name = objname)
    if(detail.save){detail$SigBCi=RC1.m}
    if(detail.null)
    {
      detail.rand$BC.obs=lapply(1:length(RC1res),function(i){RC1res[[i]]$BC.obs})
      detail.rand$BCi.rand=lapply(1:length(RC1res),function(i){RC1res[[i]]$rand})
    }
    gc()
    xx=RC1.m
  }
  if(taxa.rand.scale %in% c("across.all","both"))
  {
    if(nrow(comm)>150){big.method="loop"}else{big.method="no"}
    RC2res=iCAMP::RC.bin.bigc(com = comm,sp.bin = sp.bin,rand = rand,nworker = nworker,
                              memory.G = memory.G,big.method = big.method,
                              unit.sum = unit.sum,meta.ab=meta.ab,
                              sig.index=sig.ind.tax,detail.null=detail.null,
                              output.bray=detail.null,taxo.metric=taxo.metric,
                              transform.method=transform.method, logbase=logbase,
                              dirichlet=dirichlet)
    if(sig.ind.tax=="Confidence"){objname=paste0("C",taxo.metric,"a")}else if(sig.ind.tax=="SES"){objname=paste0("SES",taxo.metric,"a")}else if(sig.ind.tax=="RC"){objname=paste0("RC",taxo.metric,"a")}
    RC2.m=iCAMP::dist.bin.3col(RC2res$index, obj.name = objname)
    if(detail.save){detail$SigBCa=RC2.m}
    if(detail.null)
    {
      detail.rand$BC.obs=RC2res$BC.obs
      detail.rand$BCa.rand=RC2res$rand
    }
    gc()
    xx=RC2.m
  }
  
  bin.weight=data.frame(matrix(0,nrow = nrow(xx),ncol=bin.num+2))
  bin.weight[,1:2]=xx[,1:2]
  com.samp.sum=rowSums(comm)
  bin.weight[,3:(bin.num+2)]=((bin.abund[match(bin.weight[,1],rownames(bin.abund)),]/com.samp.sum[match(bin.weight[,1],rownames(comm))])
                              +(bin.abund[match(bin.weight[,2],rownames(bin.abund)),]/com.samp.sum[match(bin.weight[,2],rownames(comm))]))/2
  colnames(bin.weight)=c("samp1","samp2",paste("bin",bin.lev,sep = ""))
  if(detail.save){detail$bin.weight=bin.weight}
  
  # 8 # processes
  # 8.1 # identify governing processes
  
  
  if(sig.ind.phy=="SES"){sig.phy.cut<-sig.phy2.cut<-ses.cut}else if(sig.ind.phy=="RC"){sig.phy.cut<-sig.phy2.cut<-rc.cut}else{sig.phy.cut<-sig.phy2.cut<-conf.cut}
  if(sig.ind.tax=="SES"){sig.tax.cut<-ses.cut}else if(sig.ind.tax=="RC"){sig.tax.cut<-rc.cut}else{sig.tax.cut<-conf.cut}  
  res=integer(0)
  if((phylo.rand.scale %in% c("within.bin","both"))&(taxa.rand.scale %in% c("within.bin","both")))
  {
    if(phylo.metric %in% c("bMPD","bNRI","both"))
    {
      qp.R1R1=iCAMP::qp.bin.js(sig.phy.bin=bNRI1.m, sig.phy2.bin=NULL, sig.tax.bin=RC1.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMPDiC",taxo.name,"i")}else if(sig.index=="SES.RC"){methodname="bNRIiRCi"}else if(sig.index=="SES"){methodname=paste0("bNRIiSES",taxo.metric,"i")}else if(sig.index=="RC"){methodname=paste0("RCbMPDiRC",taxo.metric,"i")}
      res=c(res,list(qp.R1R1))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.R1R1,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric %in% c("bMNTD","bNTI","both"))
    {
      qp.T1R1=iCAMP::qp.bin.js(sig.phy.bin=bNTI1.m, sig.phy2.bin=NULL, sig.tax.bin=RC1.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDiC",taxo.name,"i")}else if(sig.index=="SES.RC"){methodname="bNTIiRCi"}else if(sig.index=="SES"){methodname=paste0("bNTIiSES",taxo.metric,"i")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDiRC",taxo.metric,"i")}
      res=c(res,list(qp.T1R1))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.T1R1,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric=="both")
    {
      qp.TR1R1=iCAMP::qp.bin.js(sig.phy.bin=bNTI1.m, sig.phy2.bin=bNRI1.m, sig.tax.bin=RC1.m,
                                bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                                sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                                check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDMPDiC",taxo.name,"i")}else if(sig.index=="SES.RC"){methodname="bNTINRIiRCi"}else if(sig.index=="SES"){methodname=paste0("bNTINRIiSES",taxo.metric,"i")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDMPDiRC",taxo.metric,"i")}
      res=c(res,list(qp.TR1R1))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.TR1R1,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
  }
  
  if((phylo.rand.scale %in% c("within.bin","both"))&(taxa.rand.scale %in% c("across.all","both")))
  {
    if(phylo.metric %in% c("bMPD","bNRI","both"))
    {
      qp.R1R2=iCAMP::qp.bin.js(sig.phy.bin=bNRI1.m, sig.phy2.bin=NULL, sig.tax.bin=RC2.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMPDiC",taxo.name,"a")}else if(sig.index=="SES.RC"){methodname="bNRIiRCa"}else if(sig.index=="SES"){methodname=paste0("bNRIiSES",taxo.metric,"a")}else if(sig.index=="RC"){methodname=paste0("RCbMPDiRC",taxo.metric,"a")}
      res=c(res,list(qp.R1R2))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.R1R2,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric %in% c("bMNTD","bNTI","both"))
    {
      qp.T1R2=iCAMP::qp.bin.js(sig.phy.bin=bNTI1.m, sig.phy2.bin=NULL, sig.tax.bin=RC2.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDiC",taxo.name,"a")}else if(sig.index=="SES.RC"){methodname="bNTIiRCa"}else if(sig.index=="SES"){methodname=paste0("bNTIiSES",taxo.metric,"a")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDiRC",taxo.metric,"a")}
      res=c(res,list(qp.T1R2))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.T1R2,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric=="both")
    {
      qp.TR1R2=iCAMP::qp.bin.js(sig.phy.bin=bNTI1.m, sig.phy2.bin=bNRI1.m, sig.tax.bin=RC2.m,
                                bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                                sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                                check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDMPDiC",taxo.name,"a")}else if(sig.index=="SES.RC"){methodname="bNTINRIiRCa"}else if(sig.index=="SES"){methodname=paste0("bNTINRIiSES",taxo.metric,"a")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDMPDiRC",taxo.metric,"a")}
      res=c(res,list(qp.TR1R2))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.TR1R2,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
  }
  
  if((phylo.rand.scale %in% c("across.all","both"))&(taxa.rand.scale %in% c("within.bin","both")))
  {
    if(phylo.metric %in% c("bMPD","bNRI","both"))
    {
      qp.R2R1=iCAMP::qp.bin.js(sig.phy.bin=bNRI2.m, sig.phy2.bin=NULL, sig.tax.bin=RC1.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMPDaC",taxo.name,"i")}else if(sig.index=="SES.RC"){methodname="bNRIaRCi"}else if(sig.index=="SES"){methodname=paste0("bNRIaSES",taxo.metric,"i")}else if(sig.index=="RC"){methodname=paste0("RCbMPDaRC",taxo.metric,"i")}
      res=c(res,list(qp.R2R1))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.R2R1,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric %in% c("bMNTD","bNTI","both"))
    {
      qp.T2R1=iCAMP::qp.bin.js(sig.phy.bin=bNTI2.m, sig.phy2.bin=NULL, sig.tax.bin=RC1.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDaC",taxo.name,"i")}else if(sig.index=="SES.RC"){methodname="bNTIaRCi"}else if(sig.index=="SES"){methodname=paste0("bNTIaSES",taxo.metric,"i")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDaRC",taxo.metric,"i")}
      res=c(res,list(qp.T2R1))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.T2R1,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric=="both")
    {
      qp.TR2R1=iCAMP::qp.bin.js(sig.phy.bin=bNTI2.m, sig.phy2.bin=bNRI2.m, sig.tax.bin=RC1.m,
                                bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                                sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                                check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDMPDaC",taxo.name,"i")}else if(sig.index=="SES.RC"){methodname="bNTINRIaRCi"}else if(sig.index=="SES"){methodname=paste0("bNTINRIaSES",taxo.metric,"i")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDMPDaRC",taxo.metric,"i")}
      res=c(res,list(qp.TR2R1))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.TR2R1,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
  }
  
  
  if((phylo.rand.scale %in% c("across.all","both"))&(taxa.rand.scale %in% c("across.all","both")))
  {
    if(phylo.metric %in% c("bMPD","bNRI","both"))
    {
      qp.R2R2=iCAMP::qp.bin.js(sig.phy.bin=bNRI2.m, sig.phy2.bin=NULL, sig.tax.bin=RC2.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMPDaC",taxo.name,"a")}else if(sig.index=="SES.RC"){methodname="bNRIaRCa"}else if(sig.index=="SES"){methodname=paste0("bNRIaSES",taxo.metric,"a")}else if(sig.index=="RC"){methodname=paste0("RCbMPDaRC",taxo.metric,"a")}
      res=c(res,list(qp.R2R2))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.R2R2,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric %in% c("bMNTD","bNTI","both"))
    {
      qp.T2R2=iCAMP::qp.bin.js(sig.phy.bin=bNTI2.m, sig.phy2.bin=NULL, sig.tax.bin=RC2.m,
                               bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                               sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                               check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDaC",taxo.name,"a")}else if(sig.index=="SES.RC"){methodname="bNTIaRCa"}else if(sig.index=="SES"){methodname=paste0("bNTIaSES",taxo.metric,"a")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDaRC",taxo.metric,"a")}
      res=c(res,list(qp.T2R2))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.T2R2,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
    if(phylo.metric=="both")
    {
      qp.TR2R2=iCAMP::qp.bin.js(sig.phy.bin=bNTI2.m, sig.phy2.bin=bNRI2.m, sig.tax.bin=RC2.m,
                                bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                                sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                                check.name=FALSE)
      if(sig.index=="Confidence"){methodname=paste0("CbMNTDMPDaC",taxo.name,"a")}else if(sig.index=="SES.RC"){methodname="bNTINRIaRCa"}else if(sig.index=="SES"){methodname=paste0("bNTINRIaSES",taxo.metric,"a")}else if(sig.index=="RC"){methodname=paste0("RCbMNTDMPDaRC",taxo.metric,"a")}
      res=c(res,list(qp.TR2R2))
      names(res)[length(res)]=methodname
      if(qp.save){utils::write.csv(qp.TR2R2,file = paste0(output.wd,"/", paste(c(prefix,"process",methodname,"csv"),collapse = ".")))}
    }
  }
  if(detail.save)
  {
    detail$processes=res
    if(is.null(unit.sum)){unit.sum.sv="asNULL"}else{unit.sum.sv=mean(unit.sum)}
    if(is.null(output.wd)){output.wdout="Not_Specified"}else{output.wdout=output.wd}
    if(is.null(transform.method)){transform.method=NA}
    detail$setting=data.frame(ds=ds, pd.cut=pd.cut, max.pd=maxpd, sp.check=sp.check,
                              phylo.rand.scale=phylo.rand.scale,taxa.rand.scale=taxa.rand.scale,
                              phylo.metric=phylo.metric,sig.index=sig.index,bin.size.limit=bin.size.limit,nworker=nworker,memory.G=memory.G,
                              rtree.save=rtree.save,detail.save=detail.save,qp.save=qp.save,detail.null=detail.null,ignore.zero=ignore.zero,
                              output.wd=output.wdout,correct.special=correct.special,
                              unit.sum.mean=unit.sum.sv,special.method=special.method,
                              ses.cut = ses.cut,rc.cut = rc.cut,conf.cut=conf.cut,
                              omit.option=omit.option[1],taxo.metric=taxo.metric, transform.method=transform.method,
                              logbase=logbase, dirichlet=dirichlet,time=format(Sys.time()-t1))
    detail$comm=comm
    res=c(res,list(detail=detail))
  }
  if(detail.null)
  {
    res=c(res,list(rand=detail.rand,special.crct=detail.special))
  }
  if(detail.save & (!is.null(output.wd))){save(res,file = paste(output.wd,"/",prefix,".iCAMP.",sig.index,".detail.rda",sep=""))}
  print(format(Sys.time()-t1))
  res
}