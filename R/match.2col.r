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
  
  if(!is.null(name.check)){name.check=as.matrix(name.check);nck2c=name.check}else{nck2c=as.matrix(check.list[[1]][,1:2])}
  ckexm=sapply(1:length(check.list),function(i){sum(nck2c[,1]!=as.vector(check.list[[i]][,1]))+sum(nck2c[,2]!=as.vector(check.list[[i]][,2]))})
  if(sum(ckexm)==0 & (!rerank)){
    if(!silent){message("All match very well. Nothing needs change.")}
    res=check.list
  }else{
    if(!silent){if(sum(ckexm)==0){message("All match very well.")}else{message(paste(names(check.list)[which(ckexm!=0)],collapse = ", ")," NOT exactly match.")}}
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
      name.int=intersect(namechk,name.int)
      if(rerank){name.int=sort(name.int);name.check=name.chk}
      name2col=name.check[match(name.int,namechk),1:2]
    }else{
      name2col=check.list[[1]][match(name.int,namel[[1]]),1:2]
    }
    
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
        if(!silent) message(names(check.list)[[i]]," has the same comparisons.")
      }
      if(ncol(check.list[[i]])>2)
      {
        res[[i]]=cbind(name2col,check.list[[i]][match(name.int,namel[[i]]),3:ncol(check.list[[i]]),drop=FALSE])
      }else{
        res[[i]]=name2col
      }
      rownames(res[[i]])=c()
    }
    names(res)=names(check.list)
  }
  res
}