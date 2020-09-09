icamp.boot<-function(icamp.result,treat,rand.time=1000,
                     compare=TRUE,silent=FALSE,
                     between.group=FALSE,ST.estimation=FALSE)
{
  # 1 # match ID
  res=icamp.result
  treat=treat[which(rownames(treat) %in% as.vector(as.matrix(res[,1:2]))),1,drop=FALSE]
  trt.lev=unique(treat[,1])
  
  # 2 # importance in each group
  procn=colnames(res[3:ncol(res)])
  if(ST.estimation){procn=c(procn,"Stochasticity")}
  id2ca=paste0(res[,1],"__",res[,2])
  id2cb=paste0(res[,2],"__",res[,1])
  
  pi.bt=lapply(1:length(trt.lev),
               function(i)
               {
                 trti=trt.lev[i]
                 if(!silent) message("bootstrapping treatment i=",i," treat=",trti,". ",date())
                 sampi=rownames(treat)[treat[,1]==trti]
                 idi=which((res[,1] %in% sampi) & (res[,2] %in% sampi))
                 pij.obs=sapply(3:ncol(res),function(j){mean(res[idi,j])})
                 if(ST.estimation){pij.obs=c(pij.obs,sum(pij.obs[3:length(pij.obs)]))}
                 pij.bt=t(sapply(1:rand.time,
                                 function(x)
                                 {
                                   sampik=sample(sampi,size = length(sampi),replace = TRUE)
                                   c2ik=utils::combn(sampik,2)
                                   c2ika=paste0(c2ik[1,],"__",c2ik[2,])
                                   ida=match(c2ika,id2ca)
                                   idb=match(c2ika,id2cb)
                                   ida[is.na(ida)]=idb[is.na(ida)]
                                   outx=sapply(3:ncol(res),function(j){mean(res[ida,j],na.rm = TRUE)})
                                   if(ST.estimation){outx=c(outx,sum(outx[3:length(outx)]))}
                                   outx
                                 }))
                 pibt=rbind(pij.obs,pij.bt)
                 rownames(pibt)=c("obs",paste0("boot",1:rand.time))
                 colnames(pibt)=procn
                 pibt
               })
  names(pi.bt)=trt.lev
  
  if(between.group & (length(trt.lev)>1))
  {
    cbn=utils::combn(x=length(trt.lev),m=2)
    
    bpi.bt=lapply(1:ncol(cbn),
                 function(u)
                 {
                   i=cbn[1,u]
                   j=cbn[2,u]
                   trti=trt.lev[i]
                   trtj=trt.lev[j]
                   if(!silent) message("bootstrapping between.treatment i=",i," j=",j,". ",date())
                   sampi=rownames(treat)[treat[,1]==trti]
                   sampj=rownames(treat)[treat[,1]==trtj]
                   idij=which(((res[,1] %in% sampi) & (res[,2] %in% sampj))|
                                ((res[,1] %in% sampj) & (res[,2] %in% sampi)))
                   pijk.obs=sapply(3:ncol(res),function(k){mean(res[idij,k])})
                   if(ST.estimation){pijk.obs=c(pijk.obs,sum(pijk.obs[3:length(pijk.obs)]))}
                   pijk.bt=t(sapply(1:rand.time,
                                   function(x)
                                   {
                                     sampik=sample(sampi,size = length(sampi),replace = TRUE)
                                     sampjk=sample(sampj,size = length(sampj),replace = TRUE)
                                     
                                     c2ijk=rbind(rep(sampik,each=length(sampjk)),
                                                 rep(sampjk,length(sampik)))
                                     c2ijka=paste0(c2ijk[1,],"__",c2ijk[2,])
                                     ida=match(c2ijka,id2ca)
                                     idb=match(c2ijka,id2cb)
                                     ida[is.na(ida)]=idb[is.na(ida)]
                                     outx=sapply(3:ncol(res),function(k){mean(res[ida,k],na.rm = TRUE)})
                                     if(ST.estimation){outx=c(outx,sum(outx[3:length(outx)]))}
                                     outx
                                   }))
                   pijbt=rbind(pijk.obs,pijk.bt)
                   rownames(pijbt)=c("obs",paste0("boot",1:rand.time))
                   colnames(pijbt)=procn
                   pijbt
                 })
    names(bpi.bt)=paste0(trt.lev[cbn[1,]],"_vs_",trt.lev[cbn[2,]])
    pi.bt=c(pi.bt,bpi.bt)
  }
  
  
  rbindlist<-function(alist)
  {
    if(is.vector(alist[[1]]))
    {
      lens=sapply(alist,length)
      alm=lapply(alist,function(v){c(v,rep(NA,max(lens)-length(v)))})
    }else{
      ncols=sapply(alist,ncol)
      alm=lapply(alist,function(v){vad=matrix(NA,nrow=nrow(v),ncol=(max(ncols)-ncol(v)));outv=cbind(as.matrix(v),vad);colnames(outv)=c();outv})
    }
    out=Reduce(rbind,alm)
    rownames(out)=c();colnames(out)=c()
    out
  }
  
  bt.sum=lapply(1:length(pi.bt),
                function(i)
                {
                  bti=pi.bt[[i]]
                  btis=lapply(1:ncol(bti),
                              function(j)
                              {
                                btij=bti[,j]
                                qtij=stats::quantile(btij,na.rm = TRUE)
                                bxij=grDevices::boxplot.stats(as.vector(btij[!is.na(btij)]))
                                c(btij[1],mean(btij,na.rm = TRUE),stats::sd(btij,na.rm = TRUE),qtij,bxij$stats,bxij$out)
                              })
                  btiss=rbindlist(btis)
                  data.frame(treat=names(pi.bt)[i],process=procn,btiss)
                })
  btsum=data.frame(rbindlist(bt.sum),stringsAsFactors = FALSE)
  colnames(btsum)=c("Group","Process","Observed","Mean","Stdev",
                    "Min","Quartile25","Median","Quartile75","Max",
                    "Lower.whisker","Lower.hinge","Median","Upper.hinge","Upper.whisker",
                    paste0("Outlier",1:(ncol(btsum)-15)))
  if(compare)
  {
    btsig<-function(x,y)
    {
      EPS=(.Machine$double.eps)
      obs.dif=x[1]-y[1]
      xm=as.vector(matrix(x[2:length(x)],nrow=length(x)-1,ncol = length(y)-1))
      ym=as.vector(matrix(y[2:length(y)],nrow=length(x)-1,ncol = length(y)-1,byrow = TRUE))
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
      c(rel.dif=rd,cohen.d=cohd$d,eff.size=cohd$magnitude,p.value=p)
    }
    
    compl=list();x=1
    for(i in 1:(length(pi.bt)-1))
    {
      pbti=pi.bt[[i]]
      for(j in (i+1):length(pi.bt))
      {
        if(!silent) message("Comparing i=",i," j=",j,". ",date())
        pbtj=pi.bt[[j]]
        p.ij=sapply(1:ncol(pbtj),
                    function(k)
                    {
                      btsig(x=pbti[,k],y=pbtj[,k])
                    })
        
        compl[[x]]=c(names(pi.bt)[i],names(pi.bt)[j],as.vector(p.ij))
        x=x+1
      }
    }
    if(length(compl)==1)
    {
      compm=data.frame(matrix(compl[[1]],nrow=1),stringsAsFactors = FALSE)
    }else{
      compm=data.frame(Reduce(rbind,compl),stringsAsFactors = FALSE)
    }
    rownames(compm)=c()
    indn=c("Relative.Diff","Cohen.d","Effect.Size","P.value")
    colnn=paste0(as.vector(matrix(procn,nrow=length(indn),ncol=length(procn),byrow = TRUE)),
                 "_",as.vector(matrix(indn,nrow=length(indn),ncol=length(procn))))
    colnames(compm)=c("Group1","Group2",colnn)
  }else{
    compm=NULL
  }
  list(summary=btsum,compare=compm,boot.detail=pi.bt)
}