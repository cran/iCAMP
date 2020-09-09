tree.droot<-function(tree,range=NA,nworker=4,output.path=FALSE)
{
  if(!is.rooted(tree)){stop("The tree is not rooted, so tree.droot stop. ", date())}
  edge<-tree$edge
  edge.len<-tree$edge.length
  tip.num=length(tree$tip)
  rootid=edge[1,1]
    
  if(is.na(range[1])){range=(1:max(tree$edge))[-rootid]}
  path=iCAMP::tree.path(tree,nworker,range)
  droot<-function(i,path){sum(path[[i]][[2]])}
  
  requireNamespace("parallel")
  c1<-parallel::makeCluster(nworker,type="PSOCK")
  message("Now computing dist to root. begin at ", date(),". Please wait...")
  dr<-parallel::parSapply(c1,1:length(path),droot,path=path)
  parallel::stopCluster(c1)
  gc()
  
  res=data.frame(node=range,distRoot=dr)
  res[nrow(res)+1,]=c(rootid,0)
  if(output.path)
  {
    output=list(droot=res,path=path)
  }else(output=res)
  output
}