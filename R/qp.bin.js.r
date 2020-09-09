qp.bin.js<-function(sig.phy.bin=NULL, sig.phy2.bin=NULL, sig.tax.bin=NULL, bin.weight,
                    sig.phy.cut=1.96, sig.phy2.cut=1.96, sig.tax.cut=0.95, check.name=FALSE)
{
  # v20200728 change bNTI.bin=NULL,bNRI.bin=NULL,RC.bin=NULL to sig.phy.bin sig.tax.bin
  if(check.name)
  {
    aslist<-function(a){if(is.null(a)){NULL}else{out=list(a);names(out)=deparse(substitute(a));out}}
    sampc=iCAMP::match.2col(check.list = c(aslist(sig.phy.bin),aslist(sig.phy2.bin),aslist(sig.tax.bin),aslist(bin.weight)))
    if(!is.null(sig.phy.bin)){sig.phy.bin=sampc$sig.phy.bin}
    if(!is.null(sig.phy2.bin)){sig.phy2.bin=sampc$sig.phy2.bin}
    if(!is.null(sig.tax.bin)){sig.tax.bin=sampc$sig.tax.bin}
    bin.weight=sampc$bin.weight
  }
  
  res=data.frame(matrix(0,nrow = nrow(bin.weight),ncol = 7))
  colnames(res)=c("sample1","sample2","Heterogeneous.Selection","Homogeneous.Selection",
                  "Dispersal.Limitation","Homogenizing.Dispersal","Drift.and.Others")
  res[,1:2]=bin.weight[,1:2]
  qpbinJS<-function(sig.phy=NULL,sig.phy2=NULL,sig.tax,weight,
                    sig.phy.cut=1.96, sig.phy2.cut=1.96, sig.tax.cut=0.95)
  {
    output=rep(0,5)
    total=sum(weight)
    if(is.null(sig.phy)){sig.phy=rep(0,length(weight))}
    if(is.null(sig.phy2)){sig.phy2=rep(0,length(weight))}
    
    output[1]=sum(weight[(sig.phy>sig.phy.cut | sig.phy2>sig.phy2.cut)])/total
    output[2]=sum(weight[(sig.phy<(-sig.phy.cut) | sig.phy2<(-sig.phy2.cut))])/total 
    output[3]=sum(weight[(abs(sig.phy)<=sig.phy.cut)&(abs(sig.phy2)<=sig.phy2.cut)&(sig.tax>sig.tax.cut)])/total
    output[4]=sum(weight[(abs(sig.phy)<=sig.phy.cut)&(abs(sig.phy2)<=sig.phy2.cut)&(sig.tax<(-sig.tax.cut))])/total
    output[5]=sum(weight[(abs(sig.phy)<=sig.phy.cut)&(abs(sig.phy2)<=sig.phy2.cut)&(abs(sig.tax)<=sig.tax.cut)])/total
    
    output
  }
  gt3n<-function(xx,i){if(is.null(xx)){out=NULL}else{out=xx[i,3:ncol(xx)]};out}
  res[,3:7]=t(sapply(1:nrow(bin.weight),
                     function(i)
                     {
                       qpbinJS(sig.phy=gt3n(xx=sig.phy.bin,i=i),
                               sig.phy2=gt3n(xx=sig.phy2.bin,i=i),
                               sig.tax=gt3n(xx=sig.tax.bin,i=i),
                               weight=bin.weight[i,3:ncol(bin.weight)],
                               sig.phy.cut=sig.phy.cut, sig.phy2.cut=sig.phy2.cut,
                               sig.tax.cut=sig.tax.cut)
                     }))
  res
}