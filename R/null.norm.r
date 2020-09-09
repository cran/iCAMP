null.norm<-function(icamp.output=NULL, rand.list=NULL, index.name="Test.Index",
                    p.norm.cut=0.05, detail.out=FALSE)
{
  if(!is.null(icamp.output))
  {
    if(is.null(icamp.output$rand)){stop("the icamp.output need to have rand (null values).")}
    rands=icamp.output$rand[grep(pattern = ".rand$",names(icamp.output$rand))]
  }else if(!is.null(rand.list)){
    rands=list(rand.list)
    names(rands)[1]=index.name
  }
  norm.test<-function(x)
  {
    if(stats::sd(x)==0){
      out=rep(NA,10)
      names(out)=c("Anderson.A","Cramer.W","Kolmogorov.D","ShapiroF.W","Shapiro.W",
                   "Anderson.P","Cramer.P","Kolmogorov.P","ShapiroF.P","Shapiro.P")
      }else{
        requireNamespace("nortest")
        adt=nortest::ad.test(x)
        cvt=nortest::cvm.test(x)
        kst=nortest::lillie.test(x)
        sft=nortest::sf.test(x)
        swt=stats::shapiro.test(x)
        out=c(Anderson=adt$statistic,Cramer=cvt$statistic,Kolmogorov=kst$statistic,ShapiroF=sft$statistic,Shapiro=swt$statistic,
              Anderson.P=adt$p.value,Cramer.P=cvt$p.value,Kolmogorov.P=kst$p.value,ShapiroF.P=sft$p.value,Shapiro.P=swt$p.value)
    }
    out
  }
  
  if(detail.out){detail.normtest=list()}
  sum.ntl=list()
  for(i in 1:length(rands))
  {
    randi=rands[[i]]
    indi=names(rands)[i]
    if(detail.out){detail.normtest[[i]]=list()}
    ntisum=list()
    for(j in 1:length(randi))
    {
      binj=paste0("bin",j)
      randij=randi[[j]]
      ntij=t(sapply(1:nrow(randij),
                    function(k)
                    {
                      randijk=as.vector(as.matrix(randij[k,3:ncol(randij)]))
                      norm.test(randijk)
                    }))
      if(detail.out){detail.normtest[[i]][[j]]=data.frame(Index=indi,BinID=binj,randij[,1:2,drop=FALSE],ntij,stringsAsFactors = FALSE)}
      ntijsum=1-(colSums(ntij[,grep(pattern = ".P$", colnames(ntij)),drop=FALSE]>=p.norm.cut,na.rm = TRUE)/nrow(ntij))
      names(ntijsum)=paste0(names(ntijsum),"nonorm")
      ntisum[[j]]=ntijsum
    }
    if(detail.out){names(detail.normtest[[i]])=paste0("bin",1:length(randi))}
    sum.ntl[[i]]=data.frame(Index=indi,BinID=paste0("bin",1:length(randi)),Reduce(rbind,ntisum),stringsAsFactors = FALSE)
    rownames(sum.ntl[[i]])=c()
  }
  
  sum.ntm=Reduce(rbind,sum.ntl)
  if(detail.out)
  {
    names(detail.normtest)=names(rands)
    output=list(summary=sum.ntm,P.value.cut=p.norm.cut,detail=detail.normtest)
  }else{
    output=list(summary=sum.ntm,P.value.cut=p.norm.cut)
  }
  output
}
