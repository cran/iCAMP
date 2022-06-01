qpen.test<-function(qpen.result,treat,rand.time=1000,
                    between.group=FALSE,out.detail=TRUE,silent=FALSE)
{
  if(inherits(qpen.result,"list")){qpm=qpen.result$result}else{qpm=qpen.result}
  head(qpm)
  
  # 1 # some handy functions
  bxv<-function(v)
  {
    qt=stats::quantile(v)
    names(qt)=c("min","q25","median","q75","max")
    bxp=grDevices::boxplot.stats(v)
    bxp5=bxp$stats
    names(bxp5)=c("lower.whisker","lower.hinge","median.boxplot","upper.hinge","upper.whisker")
    outlier=bxp$out
    if(length(outlier)>0){names(outlier)=paste0("Outlier",1:length(outlier))}
    c(mean=mean(v),sd=stats::sd(v),size=sum(!is.na(v)),qt,bxp5,outlier)
  }
  lbind<-function(lst)
  {
    lens=sapply(lst,length)
    mlen=max(lens)
    namem=names(lst[[which.max(lens)]])
    t(sapply(lst,function(v){out=c(v,rep(NA,mlen-length(v)));names(out)=namem;out}))
  }
  mbind<-function(lst)
  {
    lens=sapply(lst,ncol)
    mlen=max(lens)
    namem=colnames(lst[[which.max(lens)]])
    lnew=lapply(lst,function(v){out=cbind(v,matrix(NA,nrow = nrow(v),ncol=mlen-ncol(v)));colnames(out)=namem;out})
    out=Reduce(rbind,lnew)
    rownames(out)=c()
    out
  }
  qpn=c('Heterogeneous.Selection','Homogeneous.Selection',
        'Dispersal.Limitation','Homogenizing.Dispersal','Undominated')
  qprt<-function(qpv,...)
  {
    out=rep(0,length(qpn))
    names(out)=qpn
    qpvt=table(qpv)
    out[match(names(qpvt),names(out))]=qpvt
    out1=out/sum(out)
    out1[6]=sum(out1[1:2])
    out1[7]=sum(out1[3:5])
    names(out1)[6:7]=c('Selection','Stochasticity')
    out1
  }
  vto2<-function(v)
  {
    xid=expand.grid(1:length(v),1:length(v))
    xid=xid[xid$Var1>xid$Var2,]
    cbind(v[xid[,1]],v[xid[,2]])
  }
  ex2n<-function(x2m,x2n)
  {
    x2na=paste0(x2n[,1],"__",x2n[,2])
    x2nb=paste0(x2n[,2],"__",x2n[,1])
    xm2=paste0(x2m[,1],"__",x2m[,2])
    idxa=match(x2na,xm2)
    idxb=match(x2nb,xm2)
    idxa[is.na(idxa)]=idxb[is.na(idxa)]
    out=data.frame(x2n,x2m[idxa,3:ncol(x2m),drop=FALSE],stringsAsFactors = FALSE)
    colnames(out)=colnames(x2m)
    out=out[which(!is.na(idxa)),,drop=FALSE]
    rownames(out)=c()
    out
  }
  qpmean<-function(qpmx,...)
  {
    indmx=colMeans(qpmx[,3:6],na.rm = TRUE)
    c(indmx,qprt(qpmx[,7]))
  }
  
  btsig<-function(x,y)
  {
    EPS=(.Machine$double.eps)
    obs.dif=x[1]-y[1]
    xm=x[2:length(x)]
    ym=y[2:length(y)]
    xyd=xm-ym
    rd=obs.dif/(max(abs(c(x[1],y[1]))))
    cohd=cohend(treat = x,control = y)
    if(obs.dif>0)
    {
      p=sum(xyd<=EPS,na.rm = TRUE)/(sum(!is.na(xyd))+1)
    }else if(obs.dif<0){
      p=sum(xyd>=(-EPS),na.rm = TRUE)/(sum(!is.na(xyd))+1)
    }else{
      p=sum(abs(xyd)>EPS,na.rm=TRUE)/(sum(!is.na(xyd))+1)
    }
    c(mean.obs1=x[1],mean.obs2=y[1],rel.dif=rd,
      cohen.d=cohd$d,eff.size=cohd$magnitude,p.value=p)
  }
  
  # 2 # 
  if(out.detail){ind.out<-bt.out<-list()}
  bt.sum<-ind.sum<-comp<-list()
  
  for(i in 1:ncol(treat))
  {
    treati=treat[,i,drop=FALSE]
    trti.lev=unique(treati[,1])
    btouti<-obsmi<-sampir<-list()
    grpname<-character()
    for(j in 1:length(trti.lev))
    {
      if(!silent){message("---i=",i," j=",j,". ",date())}
      sampij=rownames(treat)[treati[,1]==trti.lev[j]]
      grpname[length(grpname)+1]=trti.lev[j]
      qpij=qpm[which((qpm[,1] %in% sampij) & (qpm[,2] %in% sampij)),,drop=FALSE]
      ind.sum[[length(ind.sum)+1]]=data.frame(GroupType=colnames(treat)[i],
                                              GroupName=trti.lev[j],
                                              Index=colnames(qpij)[3:6],
                                              lbind(lapply(3:6,function(x){bxv(qpij[,x])})),
                                              stringsAsFactors = FALSE)
      if(out.detail){ind.out[[length(ind.out)+1]]=data.frame(GroupType=colnames(treat)[i],
                                                             GroupName=trti.lev[j],qpij,stringsAsFactors = FALSE)}
      obsmi[[j]]=qpmean(qpij)
      
      sampir[[j]]=list()
      qpbtl=list()
      for(x in 1:rand.time)
      {
        sampir[[j]][[x]]=sample(sampij,length(sampij),replace = TRUE)
        qpijr=ex2n(qpij,vto2(sampir[[j]][[x]]))
        t=1
        while(nrow(qpijr)==0 & t<5)
        {
          t=t+1
          sampir[[j]][[x]]=sample(sampij,length(sampij),replace = TRUE)
          qpijr=ex2n(qpij,vto2(sampir[[j]][[x]]))
        }
        qpbtl[[x]]=qpmean(qpijr)
      }
      qpbt=Reduce(cbind,qpbtl)
      btoutij=cbind(obsmi[[j]],qpbt)
      colnames(btoutij)=c("obs",paste0("boot",1:rand.time))
      btouti[[j]]=btoutij
      
      if(out.detail){bt.out[[length(bt.out)+1]]=data.frame(GroupType=colnames(treat)[i],
                                                           GroupName=trti.lev[j],
                                                           Index=rownames(btoutij),
                                                           btoutij,stringsAsFactors = FALSE)}
      bt.sum[[length(bt.sum)+1]]=data.frame(GroupType=colnames(treat)[i],
                                            GroupName=trti.lev[j],
                                            Index=rownames(btoutij),
                                            lbind(lapply(1:nrow(btoutij),function(x){bxv(btoutij[x,])})),
                                            stringsAsFactors = FALSE)
      
    }
    
    if(between.group & (length(trti.lev)>1))
    {
      for(j in 1:(length(trti.lev)-1))
      {
        sampij=rownames(treat)[treati[,1]==trti.lev[j]]
        for(k in (j+1):length(trti.lev))
        {
          if(!silent){message("---i=",i," j=",j," k=",k,". ",date())}
          sampik=rownames(treat)[treati[,1]==trti.lev[k]]
          nameijk=paste0(trti.lev[j],"-vs-",trti.lev[k])
          grpname[length(grpname)+1]=nameijk
          qpijk=qpm[which(((qpm[,1] %in% sampij) & (qpm[,2] %in% sampik))|
                            ((qpm[,1] %in% sampik) & (qpm[,2] %in% sampij))),,drop=FALSE]
          
          ind.sum[[length(ind.sum)+1]]=data.frame(GroupType=colnames(treat)[i],
                                                  GroupName=nameijk,
                                                  Index=colnames(qpijk)[3:6],
                                                  lbind(lapply(3:6,function(x){bxv(qpijk[,x])})),
                                                  stringsAsFactors = FALSE)
          if(out.detail){ind.out[[length(ind.out)+1]]=data.frame(GroupType=colnames(treat)[i],
                                                                 GroupName=nameijk,qpijk,stringsAsFactors = FALSE)}
          obsmi[[length(obsmi)+1]]=qpmean(qpijk)
          
          
          qpbt=sapply(1:rand.time,
                      function(x)
                      {
                        idijk=expand.grid(sampir[[j]][[x]],sampir[[k]][[x]])
                        qpijkr=ex2n(qpijk,idijk)
                        qpmean(qpijkr)
                      })
          btoutijk=cbind(obsmi[[length(obsmi)]],qpbt)
          colnames(btoutijk)=c("obs",paste0("boot",1:rand.time))
          btouti[[length(btouti)+1]]=btoutijk
          
          if(out.detail){bt.out[[length(bt.out)+1]]=data.frame(GroupType=colnames(treat)[i],
                                                               GroupName=nameijk,
                                                               Index=rownames(btoutijk),
                                                               btoutijk,stringsAsFactors = FALSE)}
          bt.sum[[length(bt.sum)+1]]=data.frame(GroupType=colnames(treat)[i],
                                                GroupName=nameijk,
                                                Index=rownames(btoutijk),
                                                lbind(lapply(1:nrow(btoutijk),function(x){bxv(btoutijk[x,])})),
                                                stringsAsFactors = FALSE)
        }
      }
    }
    if(length(grpname)>1)
    {
      for(u in 1:(length(grpname)-1))
      {
        message("Now compare i=", i," u=",u," in ",(length(grpname)-1),". ",date())
        
        for(v in (u+1):length(grpname))
        {
          compuv=t(sapply(1:nrow(btouti[[u]]),function(x){btsig(btouti[[u]][x,],btouti[[v]][x,])}))
          
          comp[[length(comp)+1]]=data.frame(GroupType=colnames(treat)[i],
                                            Group1=grpname[u],Group2=grpname[v],
                                            Index=rownames(btouti[[u]]),compuv,
                                            stringsAsFactors = FALSE)
        }
      }
    }
  }
  message("Now summarize results into matrixes. ",date())
  requireNamespace("data.table")
  compare.out=data.frame(data.table::rbindlist(comp),stringsAsFactors = FALSE)
  output=list(obs.summary=mbind(ind.sum),boot.summary=mbind(bt.sum),compare=compare.out)
  if(out.detail)
  {
    message("Now output the bootstrapping details. ",date())
    output=c(output,list(group.results.detail=mbind(ind.out),boot.detail=mbind(bt.out)))
  }
  message("QPEN test is completed. ",date())
  output
}
