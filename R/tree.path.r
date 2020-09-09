tree.path<-function(tree,nworker=4,range=NA,cum=c("no","from.root","from.tip","both"))
{
  all.name=vector(length = max(tree$edge))
  all.name[-tree$edge[1,1]]=c(tree$tip.label,paste("I",1:(tree$Nnode-1),sep = ""))
  all.name[tree$edge[1,1]]="Root"
  edge<-tree$edge
  edge.len<-tree$edge.length
  if(is.na(range[1]))
  {
    range=1:length(tree$tip)
  }
  
  pathf<-function(i,edge,edge.len)
  {
    pathf<-list()
    j=1
    pathf[[1]]=edge[edge[,2]==i,1]
    pathf[[2]]=edge.len[edge[,2]==i]
    while(length(edge[edge[,2]==pathf[[1]][j],])!=0){
      pathf[[1]][j+1]=edge[edge[,2]==pathf[[1]][j],1]
      pathf[[2]][j+1]=edge.len[edge[,2]==pathf[[1]][j]]
      j=j+1
    }
    if(cum[1]=="from.tip")
    {
      pathf[[3]]=cumsum(pathf[[2]])
    }else if(cum[1]=="from.root")
    {
      pathf[[3]]=rev(cumsum(rev(pathf[[2]])))
    }else if(cum[1]=="both")
    {
      pathf[[3]]=rev(cumsum(rev(pathf[[2]])))
      pathf[[4]]=cumsum(pathf[[2]])
    }
    pathf
  }
  
  requireNamespace("parallel")
  c1<-parallel::makeCluster(nworker,type="PSOCK")
  message("Now computing path. begin at ", date(),". Please wait...")
  path<-parallel::parLapply(c1,range,pathf,edge=edge,edge.len=edge.len)
  parallel::stopCluster(c1)
  gc()
  names(path)<-all.name[range]
  path
}
