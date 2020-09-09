dist.3col<-function(dist)
{
# version 2015.5.17
dist=as.matrix(dist)
rowname=rownames(dist)
colname=colnames(dist)
rown=row(dist)
coln=col(dist)
dist.v=as.vector(stats::as.dist(dist))
rown.v=as.vector(stats::as.dist(rown))
coln.v=as.vector(stats::as.dist(coln))
res=data.frame(name1=rowname[rown.v],name2=colname[coln.v],dis=dist.v)
res
}