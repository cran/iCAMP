icamp.bins<-function(icamp.detail,treat=NULL,clas=NULL,silent=FALSE,
                     boot=FALSE,rand.time=1000,between.group=FALSE)
{
  # version 2019.12.11
  # version 2020.5.16  add between group test
  if("detail" %in% names(icamp.detail)){icamp.detail=icamp.detail$detail}
  comm=icamp.detail$comm
  samp.names=rownames(comm)
  phylo.rand.scale=as.character(icamp.detail$setting$phylo.rand.scale)
  taxa.rand.scale=as.character(icamp.detail$setting$taxa.rand.scale)
  phylo.metric=as.character(icamp.detail$setting$phylo.metric)
  bin.weight=icamp.detail$bin.weight
  
  sig.index=as.character(icamp.detail$setting$sig.index)
  ses.cut=as.numeric(icamp.detail$setting$ses.cut)
  rc.cut=as.numeric(icamp.detail$setting$rc.cut)
  conf.cut=as.numeric(icamp.detail$setting$conf.cut)
  
  if(sig.index=="SES.RC")
  {
    sig.phy.cut=ses.cut
    sig.tax.cut=rc.cut
  }else if(sig.index=="SES"){
    sig.phy.cut<-sig.tax.cut<-ses.cut
  }else if(sig.index=="RC"){
    sig.phy.cut<-sig.tax.cut<-rc.cut
  }else if(sig.index=="Confidence"){
    sig.phy.cut<-sig.tax.cut<-conf.cut
  }else{stop("wrong sig.index in icamp.bins.")}
  
  procn=c("HeS","HoS","DL","HD","DR")
  qpb<-function(pdbs,tdbs,sig.phy.cut,sig.tax.cut)
  {
    phybeta.bin = pdbs[,3:ncol(pdbs),drop=FALSE]
    taxbeta.bin = tdbs[,3:ncol(tdbs),drop=FALSE]
    out=matrix(NA,nrow=nrow(phybeta.bin),ncol=ncol(phybeta.bin))
    out[which(phybeta.bin>sig.phy.cut,arr.ind = TRUE)]="HeS"
    out[which(phybeta.bin<(-sig.phy.cut),arr.ind = TRUE)]="HoS"
    out[which((abs(phybeta.bin)<=sig.phy.cut)&(taxbeta.bin>sig.tax.cut),arr.ind = TRUE)]="DL"
    out[which((abs(phybeta.bin)<=sig.phy.cut)&(taxbeta.bin<(-sig.tax.cut)),arr.ind = TRUE)]="HD"
    out[which((abs(phybeta.bin)<=sig.phy.cut)&(abs(taxbeta.bin)<=sig.tax.cut),arr.ind = TRUE)]="DR"
    colnames(out)=substring(colnames(phybeta.bin),regexpr(".bin",colnames(phybeta.bin))+1)
    data.frame(pdbs[,1:2,drop=FALSE],out,stringsAsFactors = FALSE)
  }
  
  mis=grep("^SigbM.*D.*",names(icamp.detail))
  mjs=grep("^SigBC.*",names(icamp.detail))
  
  bwk<-ptk<-pt<-bptk<-brptk<-list();x=1
  ptuvs<-wtuvks<-list();y=1
  for(mi in mis)
  {
    pdbs=icamp.detail[[mi]]
    mi.name=names(icamp.detail)[mi]
    if(sig.index %in% c("SES.RC","SES"))
    {
      if(grepl("bMNTD",mi.name))
      {
        phy.name=paste0("bNTI",substring(mi.name,nchar(mi.name)))
      }else{
        phy.name=paste0("bNRI",substring(mi.name,nchar(mi.name)))
      }
    }else if(sig.index=="RC"){
      phy.name=sub("Sig","RC",mi.name)
    }else if(sig.index=="Confidence"){
      phy.name=sub("Sig","C",mi.name)
    }
    for(mj in mjs)
    {
      tdbs=icamp.detail[[mj]]
      mj.name=names(icamp.detail)[mj]
      if(sig.index %in% c("SES.RC","RC"))
      {
        tax.name=paste0("RCbray",substring(mj.name,nchar(mj.name)))
      }else if(sig.index=="SES"){
        tax.name=paste0("SESbray",substring(mj.name,nchar(mj.name)))
      }else if(sig.index=="Confidence"){
        tax.name=paste0("Cbray",substring(mj.name,nchar(mj.name)))
      }
      mij.name=paste0(phy.name,tax.name)
      
      wtuvk=qpb(pdbs=pdbs,tdbs=tdbs,sig.phy.cut=sig.phy.cut,sig.tax.cut=sig.tax.cut)
      if(is.null(treat)){treat=matrix("All",nrow=length(samp.names),ncol=1);rownames(treat)=samp.names;colnames(treat)="NoGroup"}
      ckid=match.2col(check.list = list(bin.weight=bin.weight,wtuvk=wtuvk),silent = TRUE)
      bin.weight=ckid$bin.weight
      wtuvk=ckid$wtuvk
      wtuvks[[y]]=data.frame(Method=mij.name,wtuvk,stringsAsFactors = FALSE)
      
      ptuv=sapply(procn,
                  function(ptx)
                  {
                    wtt=matrix(0,nrow=nrow(wtuvk),ncol = ncol(wtuvk)-2)
                    wtt[which(wtuvk[,3:ncol(wtuvk),drop=FALSE]==ptx,arr.ind = TRUE)]=1
                    rowSums(wtt*bin.weight[,3:ncol(bin.weight),drop=FALSE])/rowSums(bin.weight[,3:ncol(bin.weight),drop=FALSE])
                  })
      ptuv=data.frame(wtuvk[,1:2,drop=FALSE],ptuv,stringsAsFactors = FALSE)
      ptuvs[[y]]=data.frame(Method=mij.name,ptuv,stringsAsFactors = FALSE)
      y=y+1
      
      for(i in 1:ncol(treat))
      {
        trti.lev=unique(as.vector(treat[,i]))
        for(j in 1:length(trti.lev))
        {
          if(!silent) message("Now summarizing method=",mij.name," i=",i," j=",j,". ",date())
          sampij=rownames(treat)[which(treat[,i]==trti.lev[j])]
          bwij=bin.weight[which((bin.weight[,1] %in% sampij) & (bin.weight[,2] %in% sampij)),,drop=FALSE]
          ckij=iCAMP::match.2col(check.list = list(bwij=bwij,wtij=wtuvk,ptuvij=ptuv),silent = TRUE)
          bwij=ckij$bwij
          wtij=ckij$wtij
          ptuvij=ckij$ptuvij
          ptkij=sapply(3:ncol(bwij),function(k){sapply(procn,function(pn){sum(bwij[which(wtij[,k]==pn),k])/sum(bwij[,k])})})
          colnames(ptkij)=colnames(bwij)[3:ncol(bwij)]
          
          dpkij=lapply(1:ncol(ptkij),function(x){procn[which(ptkij[,x]==max(ptkij[,x]))]})
          dppkij=sapply(1:length(dpkij),function(x){paste(dpkij[[x]],collapse = "_")})
          dpikij=sapply(1:ncol(ptkij),function(x){max(ptkij[,x])})
          dpij=rbind(DominantProcess=dppkij,DominantProcessImportance=dpikij)
          if(boot)
          {
            if(length(unique(sampij))>2)
            {
              spair=sapply(1:nrow(bwij),function(x){paste(sort(as.vector(as.matrix(bwij[x,1:2]))),collapse = "__")})
              cbn=utils::combn(length(sampij),2)
              bt.ct=lapply(1:length(dpkij),function(x){rep(0,length(dpkij[[x]]))})
              trac=seq(from=1,to=rand.time,by=200)
              for(rt in 1:rand.time)
              {
                if(!silent){if(rt %in% trac) message("bootstrapping rt=",rt,". ",date())}
                sampijr=sample(sampij,size = length(sampij),replace = TRUE)
                tt=1
                while(length(unique(sampijr))<2 & tt<100)
                {
                  sampijr=sample(sampij,size = length(sampij),replace = TRUE)
                  tt=tt+1
                }
                spairy=sapply(1:ncol(cbn),function(x){paste(sort(sampijr[cbn[,x]]),collapse = "__")})
                spairy=spairy[which(spairy %in% spair)]
                bwijy=bwij[match(spairy,spair),,drop=FALSE]
                wtijy=wtij[match(spairy,spair),,drop=FALSE]
                ptkijy=sapply(3:ncol(bwijy),function(k){sapply(procn,function(pn){sum(bwijy[which(wtijy[,k]==pn),k])/sum(bwijy[,k])})})
                dpkijy=lapply(1:ncol(ptkijy),function(x){procn[which(ptkijy[,x]==max(ptkijy[,x]))]})
                bt.ct=lapply(1:length(dpkij),function(x){bt.ct[[x]]+(dpkij[[x]] %in% dpkijy[[x]])})
              }
              pv.bt=sapply(1:length(bt.ct),function(x){paste(round(1-(bt.ct[[x]]/rand.time),3),collapse = "_")})
            }else{
              pv.bt=rep(NA,length(dpkij))
            }
            dpij=rbind(dpij,DominantProcessPvalue=pv.bt)
          }
          
          ptkijm=data.frame(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=trti.lev[j],Index=c(procn,rownames(dpij)),rbind(ptkij,dpij),stringsAsFactors = FALSE)
          
          ptij=colMeans(ptuvij[,3:ncol(ptuvij),drop=FALSE])
          
          bptkij=sapply(3:ncol(bwij),function(k){sapply(procn,function(pn){sum(bwij[which(wtij[,k]==pn),k])/nrow(bwij)})})
          brptkij=bptkij/ptij
          colnames(bptkij)<-colnames(brptkij)<-colnames(bwij)[3:ncol(bwij)]
          bptkijm=data.frame(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=trti.lev[j],Process=procn,bptkij,stringsAsFactors = FALSE)
          brptkijm=data.frame(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=trti.lev[j],Process=procn,brptkij,stringsAsFactors = FALSE)
          
          rownames(ptkijm)<-rownames(bptkijm)<-rownames(brptkijm)<-c()
          
          bwk[[x]]=c(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=trti.lev[j],colMeans(bwij[,3:ncol(bwij)]))
          ptk[[x]]=ptkijm
          pt[[x]]=c(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=trti.lev[j],ptij)
          bptk[[x]]=bptkijm
          brptk[[x]]=brptkijm
          x=x+1
        }
        #################################################
        if(between.group)
        {
          for(j in 1:(length(trti.lev)-1))
          {
            sampij=rownames(treat)[which(treat[,i]==trti.lev[j])]
            
            for(v in (j+1):length(trti.lev))
            {
              if(!silent) message("Now summarizing method=",mij.name," i=",i," for between.group j=",j," v=",v,". ",date())
              sampiv=rownames(treat)[which(treat[,i]==trti.lev[v])]
              bwijv=bin.weight[which(((bin.weight[,1] %in% sampij) & (bin.weight[,2] %in% sampiv))|
                                       ((bin.weight[,1] %in% sampiv) & (bin.weight[,2] %in% sampij)))
                               ,,drop=FALSE]
              ckijv=iCAMP::match.2col(check.list = list(bwijv=bwijv,wtijv=wtuvk,ptuvijv=ptuv),silent = TRUE)
              bwijv=ckijv$bwijv
              wtijv=ckijv$wtijv
              ptuvijv=ckijv$ptuvijv
              ptkijv=sapply(3:ncol(bwijv),function(k){sapply(procn,function(pn){sum(bwijv[which(wtijv[,k]==pn),k])/sum(bwijv[,k])})})
              colnames(ptkijv)=colnames(bwijv)[3:ncol(bwijv)]
              
              
              dpkijv=lapply(1:ncol(ptkijv),function(x){procn[which(ptkijv[,x]==max(ptkijv[,x]))]})
              
              dppkijv=sapply(1:length(dpkijv),function(x){paste(dpkijv[[x]],collapse = "_")})
              dpikijv=sapply(1:ncol(ptkijv),function(x){max(ptkijv[,x])})
              dpijv=rbind(DominantProcess=dppkijv,DominantProcessImportance=dpikijv)
              if(boot)
              {
                spair=sapply(1:nrow(bwijv),function(x){paste(sort(as.vector(as.matrix(bwijv[x,1:2]))),collapse = "__")})
                cbn=rbind(rep(1:length(sampij),each=length(sampiv)),
                          rep(1:length(sampiv),length(sampij)))
                bt.ct=lapply(1:length(dpkijv),function(x){rep(0,length(dpkijv[[x]]))})
                trac=seq(from=1,to=rand.time,by=200)
                for(rt in 1:rand.time)
                {
                  if(!silent) {if(rt %in% trac) message("bootstrapping between.group rt=",rt,". ",date())}
                  sampijr=sample(sampij,size = length(sampij),replace = TRUE)
                  sampivr=sample(sampiv,size = length(sampiv),replace = TRUE)
                  spairy=sapply(1:ncol(cbn),function(x){paste(sort(c(sampijr[cbn[1,x]],sampivr[cbn[2,x]])),collapse = "__")})
                  spairy=spairy[which(spairy %in% spair)]
                  bwijvy=bwijv[match(spairy,spair),,drop=FALSE]
                  wtijvy=wtijv[match(spairy,spair),,drop=FALSE]
                  ptkijvy=sapply(3:ncol(bwijvy),function(k){sapply(procn,function(pn){sum(bwijvy[which(wtijvy[,k]==pn),k])/sum(bwijvy[,k])})})
                  dpkijvy=lapply(1:ncol(ptkijvy),function(x){procn[which(ptkijvy[,x]==max(ptkijvy[,x]))]})
                  bt.ct=lapply(1:length(dpkijv),function(x){bt.ct[[x]]+(dpkijv[[x]] %in% dpkijvy[[x]])})
                }
                pv.bt=sapply(1:length(bt.ct),function(x){paste(round(1-(bt.ct[[x]]/rand.time),3),collapse = "_")})
                
                dpijv=rbind(dpijv,DominantProcessPvalue=pv.bt)
              }
            
              ptkijvm=data.frame(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=paste0(trti.lev[j],"_vs_",trti.lev[v]),Index=c(procn,rownames(dpijv)),rbind(ptkijv,dpijv),stringsAsFactors = FALSE)
              
              ptijv=colMeans(ptuvijv[,3:ncol(ptuvijv),drop=FALSE])
              
              bptkijv=sapply(3:ncol(bwijv),function(k){sapply(procn,function(pn){sum(bwijv[which(wtijv[,k]==pn),k])/nrow(bwijv)})})
              brptkijv=bptkijv/ptijv
              colnames(bptkijv)<-colnames(brptkijv)<-colnames(bwijv)[3:ncol(bwijv)]
              bptkijvm=data.frame(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=paste0(trti.lev[j],"_vs_",trti.lev[v]),Process=procn,bptkijv,stringsAsFactors = FALSE)
              brptkijvm=data.frame(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=paste0(trti.lev[j],"_vs_",trti.lev[v]),Process=procn,brptkijv,stringsAsFactors = FALSE)
              rownames(ptkijvm)<-rownames(bptkijvm)<-rownames(brptkijvm)<-c()
            
              bwk[[x]]=c(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=paste0(trti.lev[j],"_vs_",trti.lev[v]),colMeans(bwijv[,3:ncol(bwijv)]))
              ptk[[x]]=ptkijvm
              pt[[x]]=c(Method=mij.name,GroupBasedOn=colnames(treat)[i],Group=paste0(trti.lev[j],"_vs_",trti.lev[v]),ptijv)
              bptk[[x]]=bptkijvm
              brptk[[x]]=brptkijvm
              x=x+1
            }
          }
        }
      }
    }
  }
  mReduce<-function(lis)
  {
    if(length(lis)==1){if(is.null(nrow(lis[[1]]))){out=data.frame(matrix(lis[[1]],nrow=1),stringsAsFactors = FALSE);colnames(out)=names(lis[[1]])}else{out=data.frame(lis[[1]],stringsAsFactors = FALSE)}}else{
      out=data.frame(Reduce(rbind,lis),stringsAsFactors = FALSE)
      rownames(out)=c()
    }
    out
  }
  output=list(Wtuvk=mReduce(wtuvks),Ptuv=mReduce(ptuvs),Ptk=mReduce(ptk),
              Pt=mReduce(pt),BPtk=mReduce(bptk),BRPtk=mReduce(brptk),Binwt=mReduce(bwk))
  
  if(!is.null(clas))
  {
    sp.bin=icamp.detail$taxabin$sp.bin[,3,drop=FALSE]
    sp.rab=colMeans(comm/rowSums(comm))
    clas=as.matrix(clas)
    spc=iCAMP::match.name(rn.list = list(clas=clas,sp.bin=sp.bin),lf.list = list(sp.rab=sp.rab),silent = TRUE)
    clas=spc$clas
    sp.bin=spc$sp.bin
    sp.rab=spc$sp.rab
    clas.bin=data.frame(Bin=paste0("Bin",sp.bin[,1]),TaxonRelativeAbundance=sp.rab,clas,stringsAsFactors = FALSE)
    bin.lev=unique(clas.bin$Bin)
    bin.lev=bin.lev[order(as.numeric(substring(bin.lev,4)),decreasing = FALSE)]
    
    bin.clas=t(sapply(1:length(bin.lev),
                      function(i)
                      {
                        idi=which(clas.bin$Bin==bin.lev[i])
                        idmaxi=idi[which.max(sp.rab[idi])]
                        topspi=rownames(clas.bin)[idmaxi]
                        topspab=max(sp.rab[idi],na.rm = TRUE)/sum(sp.rab[idi],na.rm = TRUE)
                        c(bin.lev[i],sum(sp.rab[idi]),topspi,topspab,clas[idmaxi,])
                      }))
    colnames(bin.clas)=c("Bin","BinRA","TopTaxonID","TopTaxonRAinBin",paste0("TopTaxon.",colnames(clas)))
    rownames(bin.clas)=c()
    
    bin.max=t(sapply(1:length(bin.lev),
                     function(i)
                     {
                       idi=which(clas.bin$Bin==bin.lev[i])
                       clasi=clas[idi,,drop=FALSE]
                       rabi=sp.rab[idi]
                       maxi=lapply(1:ncol(clasi),
                                   function(j)
                                   {
                                     clasij=clasi[,j]
                                     levij=unique(clasij)
                                     
                                     rabij=sapply(1:length(levij),function(k){sum(rabi[which(clasij==levij[k])])})
                                     
                                     idij.named=which((!grepl("unclass",levij))&(!grepl("Unclass",levij)))
                                     #idij.named=which(!grepl("unclass",tolower(levij)))
                                     if(length(idij.named)==0)
                                     {
                                       out=c(levij[which.max(rabij)],max(rabij)/sum(rabij))
                                     }else{
                                       out=c(levij[idij.named[which.max(rabij[idij.named])]],max(rabij[idij.named])/sum(rabij))
                                     }
                                     out
                                   })
                       unlist(maxi)
                     }))
    colnames(bin.max)=as.vector(matrix(c(paste0(colnames(clas),".maxNamed"),paste0(colnames(clas),".maxNamed.Percent")),nrow=2,byrow=TRUE))
    rownames(bin.max)=c()
    
    bin.clas=cbind(bin.clas,bin.max)
    bin.clas=data.frame(apply(bin.clas,2,as.character),stringsAsFactors = FALSE)
    idcha=c(2,4,(2*(1:ncol(clas)))+ncol(clas)+4)
    bin.clas[,idcha]=apply(bin.clas[,idcha],2,as.numeric)
    output=c(output, list(Bin.TopClass=bin.clas,Class.Bin=clas.bin))
  }
  output
}