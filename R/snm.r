snm<-function(comm,meta.com=NULL,taxon=NULL,alpha=0.05,simplify=FALSE)
{
  requireNamespace("minpack.lm")
  requireNamespace("Hmisc")
  requireNamespace("stats4")
  # Modified from the R code from Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development.ISME J. 2016 Mar;10(3):655-64
  
  comm=comm[,colSums(comm)>0]
  if(simplify){silent=TRUE}else{silent=FALSE}
  if(!is.null(meta.com))
  {
    old.sum=sum(meta.com)
    if(is.null(taxon))
    {
      spc=iCAMP::match.name(cn.list=list(comm=comm,meta.com=meta.com),silent = silent)
      comm=spc$comm
      meta.com=spc$meta.com
    }else{
      spc=iCAMP::match.name(cn.list=list(comm=comm,meta.com=meta.com),rn.list=list(taxon=taxon),silent = silent)
      comm=spc$comm
      meta.com=spc$meta.com
      taxon=spc$taxon
    }
  }
  #Calculate the number of individuals per community
  Nm <- rowSums(comm)
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(meta.com)){
    comm.ra=comm/Nm
    p = colMeans(comm.ra) # attention p.m=0 
  } else {
    p =colSums(meta.com)/old.sum
  }
  
  #Calculate the occurrence frequency of each taxa across communities
  comb=1*(comm>0)
  freq=colMeans(comb)
  
  #Combine
  C = merge(p, freq, by=0)
  C = C[order(C[,2]),]
  C.0 = data.frame(C)
  p = C.0[,2]
  freq = C.0[,3]
  names(p) = C.0[,1]
  names(freq) = C.0[,1]
  
  #Calculate the limit of detection
  N=round(mean(Nm)) # if the OTU table was rarefied to the same depth
  d = 1/N
  
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- minpack.lm::nlsLM(freq ~ stats::pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  if(!simplify)
  {
    m.ci <- try(stats::confint(m.fit, 'm', level=(1-alpha)))
    if(class(m.ci)=="try-error") m.ci=c(NA,NA)
  }
  
  
  ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
  if(!simplify)
  {
    sncm.LL <- function(m, sigma){
      R = freq - stats::pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
      R = stats::dnorm(R, 0, sigma)
      -sum(log(R))
    }
    m.mle <- try(stats4::mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p)))
    if(class(m.mle)=="try-error"){aic.fit<-bic.fit<-m.mle.coef<-m.maxll<-NA}else{
      ##Calculate Akaike's Information Criterion (AIC)
      aic.fit <- stats::AIC(m.mle, k=2)
      bic.fit <- stats::BIC(m.mle)
      m.mle.coef=m.mle@coef['m']
      m.maxll=m.mle@details$value
    }
  }
  
  ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
  freq.pred <- stats::pbeta(d, N*stats::coef(m.fit)*p, N*stats::coef(m.fit)*(1-p), lower.tail=FALSE)
  if(!simplify)
  {
    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  }
  pred.ci <- Hmisc::binconf(freq.pred*nrow(comm), nrow(comm), alpha=alpha, method="wilson", return.df=TRUE,include.x = TRUE,include.n = TRUE)
  
  ##Calculate AIC for binomial model
  if(!simplify)
  {
    bino.LL <- function(mu, sigma){
      R = freq - stats::pbinom(d, N, p, lower.tail=FALSE)
      R = stats::dnorm(R, mu, sigma)
      -sum(log(R))
    }
    bino.mle <- try(stats4::mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p)))
    if(class(bino.mle)=="try-error"){aic.bino<-bic.bino<-bino.maxll<-NA}else{
      aic.bino <- stats::AIC(bino.mle, k=2)
      bic.bino <- stats::BIC(bino.mle)
      bino.maxll=bino.mle@details$value
    }
    
    
    ##Goodness of fit for binomial model
    bino.pred <- stats::pbinom(d, N, p, lower.tail=FALSE)
    Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
    
    bino.pred.ci <- Hmisc::binconf(bino.pred*nrow(comm), nrow(comm), alpha=alpha, method="wilson", return.df=TRUE,include.x = TRUE,include.n = TRUE)
    
    ##Calculate AIC for Poisson model
    pois.LL <- function(mu, sigma){
      R = freq - stats::ppois(d, N*p, lower.tail=FALSE)
      R = stats::dnorm(R, mu, sigma)
      -sum(log(R))
    }
    pois.mle <- try(stats4::mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p)))
    if(class(pois.mle)=="try-error"){aic.pois<-bic.pois<-pois.maxll<-NA }else{
      aic.pois <- stats::AIC(pois.mle, k=2)
      bic.pois <- stats::BIC(pois.mle)
      pois.maxll<-pois.mle@details$value
    }
    
    ##Goodness of fit for Poisson model
    pois.pred <- stats::ppois(d, N*p, lower.tail=FALSE)
    Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
    
    pois.pred.ci <- Hmisc::binconf(pois.pred*nrow(comm), nrow(comm), alpha=alpha, method="wilson", return.df=TRUE,include.x = TRUE,include.n = TRUE)
    
    ##Results
    fitstats <- data.frame(m=numeric(), m.ci.2.5=numeric(),m.ci.97.5=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(stats::coef(m.fit), m.ci[1], m.ci[2], m.mle.coef, m.maxll, bino.maxll, pois.maxll, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(comm), length(p), d)
    
    A <- cbind(p, freq, freq.pred, pred.ci$Lower,pred.ci$Upper, bino.pred, bino.pred.ci$Lower,bino.pred.ci$Upper,pois.pred,pois.pred.ci$Lower,pois.pred.ci$Upper)
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr','pois.pred','pois.lwr','pois.upr')
    if(is.null(taxon)){
      B <- A[order(A[,1]),]
    } else {
      B <- merge(A, taxon, by=0, all=TRUE)
      row.names(B) <- B[,1]
      B <- B[,-1]
      B <- B[order(B[,1]),]
    }
    
    type=rep("Neutral",nrow(B))
    type[which(B$freq<B$pred.lwr)]="Below"
    type[which(B$freq>B$pred.upr)]="Above"
    B=data.frame(B,type=type)
    #save.file(B,prefix=prefix,filename = "SloanNeutralModel.eachOTU")
    type.lev=c("Neutral","Below","Above")
    type.sp=sapply(type.lev,function(tn){sum(B$type==tn)/nrow(B)})
    type.num=sapply(type.lev,function(tn){sum(B$p[which(B$type==tn)])})/sum(B$p)
  }else{
    type=rep("Neutral",length(freq))
    type[which(freq<pred.ci$Lower)]="Below"
    type[which(freq>pred.ci$Upper)]="Above"
    type.lev=c("Neutral","Below","Above")
    type.sp=sapply(type.lev,function(tn){sum(type==tn)/length(freq)})
    type.num=sapply(type.lev,function(tn){sum(p[which(type==tn)])})/sum(p)
    sp.name=lapply(type.lev,function(tn){names(freq)[which(type==tn)]})
    names(sp.name)=type.lev
  }
  if(simplify){list(type.uw=type.sp,type.wt=type.num,sp.names=sp.name)}else{list(stats=fitstats,detail=B,type.uw=type.sp,type.wt=type.num)}
}

snm.boot<-function(comm,rand=1000,meta.com=NULL,taxon=NULL,alpha=0.05,detail=TRUE)
{
  if(detail)
  {
    output=iCAMP::snm(comm=comm,meta.com = meta.com,taxon = taxon,alpha=alpha,simplify = FALSE)
  }else{
    output=NULL
  }
  obs.snm=iCAMP::snm(comm=comm,meta.com = meta.com,taxon = taxon,alpha=alpha,simplify = TRUE)
  obs=unlist(c(obs.snm$type.uw,obs.snm$type.wt))
  
  bt=data.frame(t(sapply(1:rand,
                         function(i)
                         {
                           idr=sample(1:nrow(comm),nrow(comm),replace = TRUE)
                           comr=comm[idr,,drop=FALSE]
                           out=try(iCAMP::snm(comm=comr,meta.com = meta.com,taxon = taxon,alpha=alpha,simplify = TRUE))
                           iter.n=1
                           while(class(out)=="try-error"&iter.n<100)
                           {
                             idr=sample(1:nrow(comm),nrow(comm),replace = TRUE)
                             comr=comm[idr,,drop=FALSE]
                             out=try(iCAMP::snm(comm=comr,meta.com = meta.com,taxon = taxon,alpha=alpha,simplify = TRUE))
                             iter.n=iter.n+1
                           }
                           if(class(out)=="try-error"){output=rep(NA,length(obs));warning("snm error, NAs generated")}else{
                             output=unlist(c(out$type.uw,out$type.wt))
                           }
                           output
                         })))
  
  sumv<-function(v)
  {
    qt4=stats::quantile(v)
    names(qt4)=c("min","q25","median","q75","max")
    bxp=grDevices::boxplot.stats(v)
    bxp4=c(bxp$stats)
    names(bxp4)=c("LowerWhisker","LowerHinge","Median.plot","HigherHinge","HigherWhisker")
    if(length(bxp$out)>0){bxpout=bxp$out;names(bxpout)=paste0("Outlier",1:length(bxpout));bxp4=c(bxp4,bxpout)}
    out=c(mean=mean(v,na.rm = TRUE),stdev=stats::sd(v,na.rm = TRUE),qt4,bxp4)
    out
  }
  cbindl<-function(l)
  {
    mlen=sapply(l,length)
    mn=names(l[[which(mlen==max(mlen))[1]]])
    ml=lapply(l,function(v){out=c(v,rep(NA,max(mlen)-length(v)));names(out)=mn;out})
    out=Reduce(cbind,ml)
    colnames(out)=names(l)
    out
  }
  sumbt=lapply(1:ncol(bt),function(v){sumv(bt[,v])})
  sumbtm=cbindl(sumbt)
  summ=rbind(obs,sumbtm)
  names(obs)<-colnames(bt)<-colnames(summ)<-c(paste0(names(obs.snm$type.uw),".uw"),paste0(names(obs.snm$type.wt),".wt"))
  c(output,list(summary=summ,rand=bt))
}


snm.comm<-function(comm,treat=NULL,meta.coms=NULL,meta.com=NULL,meta.group=NULL,
                   rand=1000,taxon=NULL,alpha=0.05,two.tail=TRUE,output.detail=TRUE)
{
  aslist<-function(a){if(is.null(a)){NULL}else{out=list(a);names(out)=deparse(substitute(a));out}}
  if(is.null(treat)){meta.coms<-meta.group<-NULL}
  if(!is.null(meta.coms)){meta.group<-meta.com<-NULL}
  if(!is.null(meta.com)){meta.group=NULL}
  sampc=iCAMP::match.name(rn.list = c(aslist(comm),aslist(treat),aslist(meta.group)))
  comm=sampc$comm
  if(!is.null(treat)) treat=sampc$treat
  if(!is.null(meta.group)) meta.group=sampc$meta.group
  spc=iCAMP::match.name(cn.list = c(aslist(comm),aslist(meta.com)),rn.list=aslist(taxon))
  comm=spc$comm
  if(!is.null(meta.com)) meta.com=spc$meta.com
  if(!is.null(taxon)) taxon=spc$taxon
  
  
  if(is.null(treat))
  {
    if(is.null(meta.com)){meta.com=comm}
  }else{
    if(!is.null(meta.group))
    {
      if(ncol(treat)<ncol(meta.group)){
        treat=cbind(treat,meta.group[,(ncol(treat)+1):ncol(meta.group),drop=FALSE])
        warning("Since meta.group has more columns than treat, additional columns of meta.group were used as treatment.")
      }else if(ncol(treat)>ncol(meta.group)){
        meta.group=cbind(meta.group,treat[,(ncol(meta.group)+1):ncol(treat),drop=FALSE])
        warning("Since meta.group has less columns than treat, additional columns of treat were used as meta.group.")
      }
    }
    if(is.null(meta.coms)){meta.coms=list()}else{if(length(meta.coms)!=ncol(treat)){stop("meta.coms does not match treat.")}}
    for(i in 1:ncol(treat))
    {
      trti.lev=unique(treat[,i])
      if(length(meta.coms)>=i){if(length(meta.coms[[i]])!=length(trti.lev)){stop("meta.coms does not match treat.")}}else{
        meta.coms[[i]]=list()
        for(j in 1:length(trti.lev))
        {
          if(is.null(meta.group)){if(is.null(meta.com)){meta.coms[[i]][[j]]=comm}else{meta.coms[[i]][[j]]=meta.com}}else{
            metaij=unique(meta.group[which(treat[,i]==trti.lev[j]),i])
            if(length(metaij)>1){stop("meta.group setting is wrong. The same treatment should be in the same metacommunity.")}
            sampij=rownames(meta.group)[which(meta.group[,i]==metaij)]
            if(!is.null(meta.com))
            {
              if(sum(!(sampij %in% rownames(meta.com)))>0){warning("meta.com rownames did not match meta.group.")}
              meta.coms[[i]][[j]]=meta.com[which(rownames(meta.com) %in% sampij),,drop=FALSE]
            }else{
              meta.coms[[i]][[j]]=comm[which(rownames(comm) %in% sampij),,drop=FALSE]
            }
          }
        }
      }
    }
  }
  
  statses=list()
  plot.detail=list()
  ratio.sum=list()
  rands=list()
  k=1;m=1
  pvalue=list()
  
  mx=matrix(NA,nrow=rand,ncol=rand)
  idx=as.vector(col(mx))
  idy=as.vector(row(mx))
  
  if(is.null(treat))
  {
    snmi=iCAMP::snm.boot(comm,rand=rand,meta.com = meta.com,taxon=taxon,alpha = alpha,detail = TRUE)
    stats.out=data.frame(treat.type="NA",treatment.id="All",snmi$stats)
    detaili=snmi$detail
    plot.detail.out=data.frame(OTU=rownames(detaili),treat.type=rep("NA",nrow(detaili)),treatment.id=rep("All",nrow(detaili)),detaili)
    ratio.sum.out=data.frame(index=rownames(snmi$summary),treat.type=rep("NA",nrow(snmi$summary)),treatment.id=rep("All",nrow(snmi$summary)),snmi$summary)
    pvalue.out=NULL
    rownames(plot.detail.out)<-rownames(ratio.sum.out)<-c()
  }else{
    for(i in 1:ncol(treat))
    {
      trt.lev=unique(treat[,i])
      rands[[i]]=list()
      for(j in 1:length(trt.lev))
      {
        trtij=trt.lev[j]
        sampij=rownames(treat)[treat[,i]==trtij]
        comij=comm[which(rownames(comm) %in% sampij),,drop=FALSE]
        snmij=iCAMP::snm.boot(comm = comij,meta.com=meta.coms[[i]][[j]],taxon=taxon,alpha=alpha, detail = TRUE,rand = rand)
        statses[[k]]=cbind(colnames(treat)[i],trtij,snmij$stats)
        detailij=snmij$detail
        plot.detail[[k]]=data.frame(OTU=rownames(detailij),treat.type=rep(colnames(treat)[i],nrow(detailij)),treatment.id=trtij,detailij)
        rownames(plot.detail[[k]])=c()
        ratio.sum[[k]]=data.frame(index=rownames(snmij$summary),treat.type=rep(colnames(treat)[i],nrow(snmij$summary)),treatment.id=rep(trtij,nrow(snmij$summary)),snmij$summary)
        rownames(ratio.sum[[k]])=c()
        rands[[i]][[j]]=snmij$rand
        k=k+1
      }
      names(rands[[i]])=trt.lev
      
      for(x in 1:(length(trt.lev)-1))
      {
        for(y in (x+1):length(trt.lev))
        {
          randx=rands[[i]][[x]][idx,,drop=FALSE]
          randy=rands[[i]][[y]][idy,,drop=FALSE]
          p=((colSums(randx>randy,na.rm = TRUE)+0.5*colSums(randx==randy,na.rm = TRUE))/colSums(!is.na(randx)))
          p[p>0.5]=1-p[p>0.5]
          if(two.tail){p=p*2}
          pvalue[[m]]=c(colnames(treat)[i],trt.lev[x],trt.lev[y],p)
          m=m+1
        }
      }
    }
    names(rands)=colnames(treat)
    stats.out=data.frame(Reduce(rbind,statses),stringsAsFactors = FALSE)
    colnames(stats.out)[1:2]=c("treat.type","treatment.id")
    plot.detail.out=Reduce(rbind,plot.detail)
    ratio.sum.out=Reduce(rbind,ratio.sum)
    pvalue.out=data.frame(matrix(Reduce(rbind,pvalue),nrow=length(pvalue)),stringsAsFactors = FALSE)
    rownames(pvalue.out)=c()
    colnames(pvalue.out)=c("treat.type","treatment1","treatment2",names(pvalue[[1]])[-(1:3)])
  }
  
  output=list(stats=stats.out,plot.detail=plot.detail.out,ratio.summary=ratio.sum.out,pvalues=pvalue.out)
  if(output.detail){output=c(output,list(boot.detail=rands))}
  output
}


