icamp.cate<-function(icamp.bins.result,comm,cate,treat=NULL,silent=FALSE,
                     between.group=FALSE)
{
  Wtuvk=icamp.bins.result$Wtuvk
  taxbin=icamp.bins.result$Class.Bin[,1,drop=FALSE]
  comra=comm/rowSums(comm)
  procn=c("HeS","HoS","DL","HD","DR")
  
  # function idnx is to make this function compatible to most of previous versions.
  idnx<-function(m,x=1){m[,which(tolower(colnames(m)) %in% paste0(c("name","samp","sample"),x))[1]]}
  
  cate.lev=unique(cate[,1])
  BPtuvxl<-Ptuvxl<-list()
  for(i in 1:length(cate.lev))
  {
    if(!silent){message("Calculating for each category in each turnover, i=",i," in ",length(cate.lev),". ",date())}
    spi=rownames(cate)[cate[,1]==cate.lev[i]]
    bini=tolower(taxbin[match(spi,rownames(taxbin)),1])
    idi=match(bini,tolower(colnames(Wtuvk)))
    wtuvki=Wtuvk[,idi,drop=FALSE]
    rabi=(comra[match(idnx(Wtuvk,1),rownames(comra)),match(spi,colnames(comra)),drop=FALSE]+
            comra[match(idnx(Wtuvk,2),rownames(comra)),match(spi,colnames(comra)),drop=FALSE])/2
    
    BPtuvxi=sapply(1:length(procn),function(j){rowSums((wtuvki==procn[j])*rabi,na.rm = TRUE)})
    BPtuvxi[which(rowSums(rabi)==0),]=NA
    Ptuvxi=BPtuvxi/rowSums(rabi)
    Ptuvxi[which(rowSums(rabi)==0),]=NA
    BPtuvxl[[i]]=cbind(BPtuvxi,rowSums(BPtuvxi[,3:5]))
    Ptuvxl[[i]]=cbind(Ptuvxi,rowSums(Ptuvxi[,3:5]))
    colnames(BPtuvxl[[i]])<-colnames(Ptuvxl[[i]])<-paste0(cate.lev[i],".",c(procn,"Stochasticity"))
  }
  BPtuvx=Reduce(cbind,BPtuvxl)
  Ptuvx=Reduce(cbind,Ptuvxl)
  BPtuvct=data.frame(Wtuvk[,1:3,drop=FALSE],BPtuvx,stringsAsFactors = FALSE)
  Ptuvct=data.frame(Wtuvk[,1:3,drop=FALSE],Ptuvx,stringsAsFactors = FALSE)
  
  if(is.null(treat))
  {
    treat=data.frame(Group=rep('All',nrow(comm)),stringsAsFactors = FALSE)
    rownames(treat)=rownames(comm)
  }
  
  method.lev=unique(Ptuvct$Method)
  BPtxm<-Ptxm<-list()
  for(m in 1:length(method.lev))
  {
    BPtxl<-Ptxl<-list()
    for(i in 1:ncol(treat))
    {
      treati=treat[,i,drop=FALSE]
      trt.lev=unique(treati[,1])
      
      Ptxi=t(sapply(1:length(trt.lev),
                    function(j)
                    {
                      sampj=rownames(treati)[treati[,1]==trt.lev[j]]
                      colMeans(Ptuvx[which((idnx(Ptuvct,1) %in% sampj)&(idnx(Ptuvct,2) %in% sampj)),,drop=FALSE],na.rm = TRUE)
                    }))
      BPtxi=t(sapply(1:length(trt.lev),
                    function(j)
                    {
                      sampj=rownames(treati)[treati[,1]==trt.lev[j]]
                      colMeans(BPtuvx[which((idnx(BPtuvct,1) %in% sampj)&(idnx(BPtuvct,2) %in% sampj)),,drop=FALSE],na.rm = TRUE)
                    }))
      namexi=trt.lev
      if(length(trt.lev)>1 & between.group)
      {
        bptxijkl<-ptxijkl<-list()
        namejkl=list()
        x=1
        for(j in 1:(length(trt.lev)-1))
        {
          sampj=rownames(treati)[treati[,1]==trt.lev[j]]
          for(k in (j+1):length(trt.lev))
          {
            sampk=rownames(treati)[treati[,1]==trt.lev[k]]
            ptxijkl[[x]]=colMeans(Ptuvx[which(((idnx(Ptuvct,1) %in% sampj)&(idnx(Ptuvct,2) %in% sampk))|
                                                ((idnx(Ptuvct,1) %in% sampk)&(idnx(Ptuvct,2) %in% sampj))),,drop=FALSE],na.rm = TRUE)
            bptxijkl[[x]]=colMeans(BPtuvx[which(((idnx(BPtuvct,1) %in% sampj)&(idnx(BPtuvct,2) %in% sampk))|
                                                ((idnx(BPtuvct,1) %in% sampk)&(idnx(BPtuvct,2) %in% sampj))),,drop=FALSE],na.rm = TRUE)
            
            namejkl[[x]]=paste0(trt.lev[j],"_vs_",trt.lev[k])
            x=x+1
          }
        }
        BPtxi=rbind(BPtxi,Reduce(rbind,bptxijkl))
        Ptxi=rbind(Ptxi,Reduce(rbind,ptxijkl))
        rownames(Ptxi)<-rownames(BPtxi)<-c()
        namexi=c(namexi,unlist(namejkl))
      }
      BPtxl[[i]]=data.frame(Method=method.lev[m],GroupBasedOn=colnames(treat)[i],Group=namexi,BPtxi,stringsAsFactors = FALSE)
      Ptxl[[i]]=data.frame(Method=method.lev[m],GroupBasedOn=colnames(treat)[i],Group=namexi,Ptxi,stringsAsFactors = FALSE)
    }
    BPtxm[[m]]=Reduce(rbind,BPtxl)
    Ptxm[[m]]=Reduce(rbind,Ptxl)
  }
  BPtx=Reduce(rbind,BPtxm)
  Ptx=Reduce(rbind,Ptxm)
  rownames(Ptx)=c()
  
  colnames(Ptuvct)[2:3]<-colnames(BPtuvct)[2:3]<-c('samp1','samp2') # to make all outputs with consistent colnames.
  
  list(Ptuvx=Ptuvct,Ptx=Ptx,Ptuvx.wt=BPtuvct,Ptx.wt=BPtx)
}