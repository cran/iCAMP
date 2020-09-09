match.2col<-function(check.list,name.check=NULL,rerank=TRUE,silent=FALSE)
{
  if(length(names(check.list))==0)
  {
    names(check.list)=paste("obj",1:length(check.list),sep = ".")
  }else{
    for(i in 1:length(check.list))
    {
      if(nchar(names(check.list)[i])==0)
      {
        names(check.list)[i]=paste("obj",i,sep = ".")
      }
    }
  }
  
  namel=list()
  for(i in 1:length(check.list))
  {
    temp=as.matrix(check.list[[i]])
    oldname=temp[,1:2,drop=FALSE]
    oldname=t(sapply(1:nrow(oldname),function(j){oldname[j,1:2][order(oldname[j,1:2])]}))
    namel[[i]]=paste(oldname[,1],oldname[,2],sep = "_")
    if(rerank)
    {
      if(!silent) message(names(check.list)[[i]]," has been re-ranked by names.")
      check.list[[i]][,1:2]=as.matrix(oldname)
      reid=order(namel[[i]])
      check.list[[i]]=check.list[[i]][reid,]
      rownames(check.list[[i]])<-c()
      namel[[i]]=namel[[i]][reid]
    }
    
    if(i==1)
    {
      name.int=namel[[i]]
    }else{
      name.int=intersect(name.int,namel[[i]])
    }
  }
  
  if(!is.null(name.check))
  {
    name.chk=t(sapply(1:nrow(name.check), function(j){name.check[j,1:2][order(name.check[j,1:2])]}))
    namechk=paste(name.chk[,1],name.chk[,2],sep = "_")
  }else{
    namechk=unique(unlist(namel))
  }
  
  name.int=intersect(name.int,namechk)
  
  res=list()
  
  for(i in 1:length(check.list))
  {
    los.id=(!(name.int %in% namel[[i]]))
    add.id=(!(namel[[i]] %in% name.int))
    mis=sum(los.id)+sum(add.id)
    if(mis>0){
      if(!silent) message("Mismatch warning: ",names(check.list)[[i]]," has ",mis," mismatched names.")
      if(!silent) message("Names expected but not appear: ",paste(name.int[los.id],collapse = ", "),".")
      if(!silent) message("Names unexpected: ", paste(namel[[i]][add.id],collapse = ", "),".")
    }else{
      if(!silent) message(names(check.list)[[i]]," matches well.")
    }
    res[[i]]=check.list[[i]][(namel[[i]] %in% name.int),]
  }
  names(res)=names(check.list)
  res
}