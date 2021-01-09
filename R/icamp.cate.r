icamp.cate<-function(icamp.bins.result,comm,cate,treat=NULL,silent=FALSE,
                     between.group=FALSE)
{
  Wtuvk=icamp.bins.result$Wtuvk
  taxbin=icamp.bins.result$Class.Bin[,1,drop=FALSE]
  comra=comm/rowSums(comm)
  procn=c("HeS","HoS","DL","HD","DR")
  
  cate.lev=unique(cate[,1])
  Ptuvxl=list()
  for(i in 1:length(cate.lev))
  {
    if(!silent){message("Calculating for each category in each turnover, i=",i," in ",length(cate.lev),". ",date())}
    spi=rownames(cate)[cate[,1]==cate.lev[i]]
    bini=tolower(taxbin[match(spi,rownames(taxbin)),1])
    idi=match(bini,tolower(colnames(Wtuvk)))
    wtuvki=Wtuvk[,idi,drop=FALSE]
    rabi=(comra[match(Wtuvk$name1,rownames(comra)),match(spi,colnames(comra)),drop=FALSE]+
            comra[match(Wtuvk$name2,rownames(comra)),match(spi,colnames(comra)),drop=FALSE])/2
    
    Ptuvxi=sapply(1:length(procn),function(j){rowSums((wtuvki==procn[j])*rabi,na.rm = TRUE)})/rowSums(rabi)
    Ptuvxi[which(rowSums(rabi)==0),]=0
    Ptuvxl[[i]]=cbind(Ptuvxi,rowSums(Ptuvxi[,3:5]))
    colnames(Ptuvxl[[i]])=paste0(cate.lev[i],".",c(procn,"Stochasticity"))
  }
  Ptuvx=Reduce(cbind,Ptuvxl)
  Ptuvct=data.frame(Wtuvk[,1:3,drop=FALSE],Ptuvx,stringsAsFactors = FALSE)
  
  if(!is.null(treat))
  {
    
    method.lev=unique(Ptuvct$Method)
    Ptxm=list()
    for(m in 1:length(method.lev))
    {
      Ptxl=list()
      for(i in 1:ncol(treat))
      {
        treati=treat[,i,drop=FALSE]
        trt.lev=unique(treati[,1])
        
        Ptxi=t(sapply(1:length(trt.lev),
                      function(j)
                      {
                        sampj=rownames(treati)[treati[,1]==trt.lev[j]]
                        colMeans(Ptuvx[which((Ptuvct$name1 %in% sampj)&(Ptuvct$name2 %in% sampj)),,drop=FALSE])
                      }))
        namexi=trt.lev
        if(length(trt.lev)>1 & between.group)
        {
          ptxijkl=list()
          namejkl=list()
          x=1
          for(j in 1:(length(trt.lev)-1))
          {
            sampj=rownames(treati)[treati[,1]==trt.lev[j]]
            for(k in (j+1):length(trt.lev))
            {
              sampk=rownames(treati)[treati[,1]==trt.lev[k]]
              ptxijkl[[x]]=colMeans(Ptuvx[which(((Ptuvct$name1 %in% sampj)&(Ptuvct$name2 %in% sampk))|
                                                  ((Ptuvct$name1 %in% sampk)&(Ptuvct$name2 %in% sampj))),,drop=FALSE])
              namejkl[[x]]=paste0(trt.lev[j],"_vs_",trt.lev[k])
              x=x+1
            }
          }
          Ptxi=rbind(Ptxi,Reduce(rbind,ptxijkl))
          rownames(Ptxi)=c()
          namexi=c(namexi,unlist(namejkl))
        }
        
        Ptxl[[i]]=data.frame(Method=method.lev[m],GroupBasedOn=colnames(treat)[i],Group=namexi,Ptxi,stringsAsFactors = FALSE)
      }
      Ptxm[[m]]=Reduce(rbind,Ptxl)
    }
    Ptx=Reduce(rbind,Ptxm)
    rownames(Ptx)=c()
  }
  list(Ptuvx=Ptuvct,Ptx=Ptx)
}