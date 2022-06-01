maxbigm<-function(m.desc,m.wd,nworker=1,rm.na=TRUE,size.limit=10000*10000)
{
  requireNamespace("bigmemory")
  
  mx=bigmemory::attach.big.matrix(dget(paste0(m.wd,"/",m.desc)))
  coln=ncol(mx)
  rown=nrow(mx)
  inv=max(floor(size.limit/rown),1)
  num=ceiling(coln/inv)
  if(num==1){ser=matrix(c(1,coln),nrow=1)}else{
    ser=cbind((0:(num-1))*inv+1,c((1:(num-1))*inv,coln))
  }
  
  findmax<-function(i,m.desc,ser,rm.na,m.wd)
  {
    requireNamespace("bigmemory")
    mx=bigmemory::attach.big.matrix(dget(paste0(m.wd,"/",m.desc)))
    mi=mx[,ser[i,1]:ser[i,2],drop=FALSE]
    gc()
    maxi=max(mi,na.rm = rm.na)
    id=which(mi==maxi,arr.ind = TRUE)
    id[,2]=id[,2]+ser[i,1]-1
    gc()
    list(maxi,id)
  }
  
  if(nworker==1)
  {
    maxs=lapply(1:num, findmax,m.desc=m.desc,ser=ser,rm.na=rm.na,m.wd=m.wd)
  }else{
    requireNamespace("parallel")
    c1<-try(parallel::makeCluster(nworker,type="PSOCK"))
    if(inherits(c1,"try-error")){c1 <- try(parallel::makeCluster(nworker, setup_timeout = 0.5))}
    if(inherits(c1,"try-error")){c1 <- parallel::makeCluster(nworker, setup_strategy = "sequential")}
    maxs<-parallel::parLapply(c1,1:num,findmax,m.desc=m.desc,ser=ser,rm.na=rm.na,m.wd=m.wd)
    parallel::stopCluster(c1)
    gc()
  }
  
  max.v=sapply(1:length(maxs), function(j){maxs[[j]][[1]]})
  maxvalue=max(max.v)
  maxid=which(max.v==maxvalue)
  maxrc=lapply(maxid, function(u){maxs[[u]][[2]]})
  maxrc=Reduce(rbind,maxrc)
  list(max.value=maxvalue,row.col=maxrc)
}