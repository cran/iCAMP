ps.bin<-function(sp.bin,sp.ra,spname.use=NULL,
                 pd.desc = "pd.desc",pd.spname,pd.wd,
                 nd.list,nd.spname=NULL,ndbig.wd=NULL,
                 cor.method=c("pearson","spearman"),
                 r.cut=0.01,p.cut=0.2,min.spn=6)
{
  requireNamespace("bigmemory")
  if(is.null(spname.use)){spname.use=intersect(pd.spname,nd.spname)}
  
  bin.lev=unique(sp.bin[,1])
  env.names=names(nd.list)
  sp.ra=as.matrix(sp.ra)[,1]
  
  if(!is.null(ndbig.wd))
  {
    ndm=list()
    for(j in 1:length(env.names))
    {
      ndm[[j]]=bigmemory::attach.big.matrix(dget(paste0(ndbig.wd,"/",nd.list[[j]])))
    }
    names(ndm)=env.names
  }else{
    ndm=nd.list
    if(is.null(nd.spname)){nd.spname=rownames(ndm[[1]])}
  }
  
  pdm=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
  
  mtm=t(sapply(1:length(bin.lev),
               function(i)
               {
                 bini=bin.lev[i]
                 spi=rownames(sp.bin)[sp.bin[,1]==bin.lev[i]]
                 spiu=spi[which(spi %in% spname.use)]
                 binabi=sum(sp.ra[which(names(sp.ra) %in% spiu)])
                 if(length(spiu)<min.spn)
                 {
                   out=c(rep(NA,(length(env.names)*2*length(cor.method))),abu=NA)
                 }else{
                   pdi=pdm[match(spiu,pd.spname),match(spiu,pd.spname)]
                   rpi=sapply(1:length(env.names),
                              function(j)
                              {
                                if(!is.null(ndbig.wd))
                                {
                                  ndij=ndm[[j]][match(spiu,nd.spname),match(spiu,nd.spname)]
                                }else{
                                  ndij=ndm[[j]][match(spiu,rownames(ndm[[j]])),match(spiu,colnames(ndm[[j]]))]
                                }
                                outijm=sapply(1:length(cor.method),
                                              function(m)
                                              {
                                                mtijm=vegan::mantel(xdis = stats::as.dist(ndij),ydis = stats::as.dist(pdi),method = cor.method[m])
                                                c(mtijm$statistic,mtijm$signif)
                                              })
                                as.vector(t(outijm))
                              })
                   out=c(as.vector(t(rpi)),abu=binabi)
                 }
                 out
               }))
  rownames(mtm)=paste0("bin",bin.lev)
  colnames(mtm)[1:(ncol(mtm)-1)]=paste0(rep(env.names,2*length(cor.method)),".",rep(rep(cor.method,each=length(env.names)),2),rep(c(".r",".p"),each=length(env.names)*length(cor.method)))
  mtm.r=mtm[,1:((ncol(mtm)-1)/2),drop=FALSE]
  mtm.p=mtm[,((ncol(mtm)-1)/2)+(1:((ncol(mtm)-1)/2)),drop=FALSE]
  bin.abu=mtm[,ncol(mtm)]
  if(length(r.cut)<length(p.cut)){r.cut=c(r.cut,rep(r.cut[length(r.cut)],length(p.cut)-length(r.cut)))}
  if(length(r.cut)>length(p.cut)){p.cut=c(p.cut,rep(p.cut[length(p.cut)],length(r.cut)-length(p.cut)))}
  RAMR=matrix(NA,nrow=((3*length(r.cut))+1),ncol = ncol(mtm.r))
  for(i in 1:length(r.cut))
  {
    sigm=(mtm.r>=r.cut[i] & mtm.p<=p.cut[i])
    RAMR[i,]=colSums(matrix(bin.abu,nrow=nrow(mtm),ncol=ncol(mtm.r))*sigm,na.rm = TRUE)
    RAMR[length(r.cut)+i,]=RAMR[i,]/sum(bin.abu,na.rm = TRUE)
    RAMR[2*length(r.cut)+i,]=colSums(mtm.r*sigm*bin.abu,na.rm = TRUE)/colSums(sigm*bin.abu,na.rm = TRUE)
  }
  RAMR[nrow(RAMR),]=colSums(mtm.r*bin.abu,na.rm = TRUE)/sum(bin.abu,na.rm = TRUE)
  colnames(RAMR)=gsub("\\.r$","",colnames(mtm.r))
  RAMRm=data.frame(index=c(rep(c("RAsig","RAsig.adj","MeanR.sig"),each=length(r.cut)),"MeanR"),
                   r.cutoff=c(rep(r.cut,3),NA),p.cutoff=c(rep(p.cut,3),NA),RAMR,stringsAsFactors = FALSE)
  list(Index=RAMRm,detail=mtm)
}