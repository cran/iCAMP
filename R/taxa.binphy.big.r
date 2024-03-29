taxa.binphy.big<-function(tree,pd.desc, pd.spname, pd.wd, 
                          outgroup.tip=NA,outgroup.rm=TRUE,
                          d.cut=NULL,ds=0.2,bin.size.limit=24,
                          nworker=4,d.cut.method=c("maxpd","maxdroot"))
{
  requireNamespace("ape")
  requireNamespace("bigmemory")
  
  pd=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
  
  if(!is.na(outgroup.tip))
  {
    tree=ape::root(tree,outgroup = outgroup.tip,r=TRUE)
    if(outgroup.rm){tree=ape::drop.tip(tree,outgroup.tip)}
  }
  if(!is.rooted(tree)){stop("Failed to root the tree. Stop. ",date())}
  
  path=iCAMP::tree.path(tree,nworker,cum="from.root")
  
  if(is.null(d.cut))
  {
    if(is.null(ds))
    {
      stop("d.cut or ds, at least one should be given.")
    }else if(d.cut.method[1]=="maxpd"){
      pdmax=iCAMP::maxbigm(m.desc = pd.desc, m.wd = pd.wd, nworker = nworker,
                    rm.na = TRUE, size.limit = 10000 * 10000)$max.value
      d.cut=(pdmax-ds)/2
    }else if(d.cut.method[1]=="maxdroot"){
      maxdroot=max(sapply(1:length(path),function(i){max(path[[i]][[3]])}))
      d.cut=maxdroot-(ds/2)
    }
  }
  
  node=sapply(1:length(path),
              function(i)
              {
                id=which(path[[i]][[3]]>=d.cut)
                if(length(id)<=1){i}else{path[[i]][[1]][max(id)-1]}
              })
  node.table=as.matrix(table(node))
  
  sp.num=length(tree$tip.label)
  sp.name=tree$tip.label
  sp.bin<-data.frame(matrix(NA,nrow=sp.num,ncol = 3))
  colnames(sp.bin)=c("bin.id.strict","bin.id.united","bin.id.new")
  rownames(sp.bin)=sp.name
  bin.sp.id<-newbin.sp.id<-list()
  sp.core<-data.frame(matrix(nrow=nrow(node.table),ncol=5))
  colnames(sp.core)=c("bin.strict.id","bin.strict.taxa.num","bin.pd.max","bin.pd.mean","bin.pd.sd")
  
  sp.bin[,1]=node
  sp.core[,1]=as.integer(rownames(node.table)) #revised 20210626
  sp.core[,2]=node.table[,1] #revised 20210626
  
  bin.sp.id=lapply(1:nrow(sp.core), function(i){rownames(sp.bin)[which(sp.bin[,1]==sp.core[i,1])]})
  names(bin.sp.id)=sp.core[,1]
  pd.state=sapply(1:nrow(sp.core),
                  function(i)
                  {
                    if(length(bin.sp.id[[i]])==1){out=c(0,0,NA)}else{
                      pdid=match(bin.sp.id[[i]],pd.spname)
                      pd.temp=stats::as.dist(pd[pdid,pdid])
                      out=c(max(pd.temp),mean(pd.temp),stats::sd(pd.temp))
                    }
                    out
                  })
  sp.core[,3:5]=t(matrix(pd.state,nrow = 3,ncol=nrow(sp.core)))
  
  # move very small bins to the nearest good bin.
  bin.good=which(sp.core[,2]>=bin.size.limit)
  bin.small=which(sp.core[,2]<bin.size.limit)

  bin.good.node=sp.core[,1][bin.good]
  bin.small.node=sp.core[,1][bin.small]
  sp.bin.temp=sp.bin[,1]
  
  if(length(bin.small)>0)
  {
    bg.path<-bg.path.len<-bg.path.cum<-bs.path<-bs.path.len<-bs.path.cum<-list()
    if(length(bin.good)==0)
    {
      bg.path=integer(0)
      bg.path.len=integer(0)
      bg.path.cum=integer(0)
    }else{
      for(u in 1:length(bin.good))
      {
        id.tip=which(sp.bin.temp==bin.good.node[u])[1]
        range=(which(path[[id.tip]][[1]]==bin.good.node[u]):length(path[[id.tip]][[1]]))
        bg.path[[u]]=path[[id.tip]][[1]][range]
        bg.path.len[[u]]=c(0,path[[id.tip]][[2]][range[-1]])
        bg.path.cum[[u]]=cumsum(bg.path.len[[u]])
      }
    }
    
    for(u in 1:length(bin.small))
    {
      id.tip=which(sp.bin.temp==bin.small.node[u])[1]
      if(id.tip<=sp.num)
      {
        bs.path[[u]]=c(id.tip,path[[id.tip]][[1]])
        bs.path.len[[u]]=c(0,path[[id.tip]][[2]])
      }else{
        range=(which(path[[id.tip]][[1]]==bin.small.node[u]):length(path[[id.tip]][[1]]))
        bs.path[[u]]=path[[id.tip]][[1]][range]
        bs.path.len[[u]]=c(0,path[[id.tip]][[2]][range[-1]])
      }
      bs.path.cum[[u]]=cumsum(bs.path.len[[u]])
    }
           
    rootid=tree$edge[1,1]
    nodes=((length(tree$tip.label)+1):max(tree$edge));nodes=nodes[which(nodes!=rootid)]
    droot=tree.droot(tree,range = nodes,nworker,output.path = TRUE)
    node.path=droot$path
    names(node.path)=nodes
    droot=droot$droot
    droot=droot[order(droot[,2],decreasing = TRUE),]
    drank.id=match(c(bin.good.node,bin.small.node),droot[,1])
    if(sum(!is.na(drank.id))==0){drank=as.integer(droot[,1])}else{drank=as.integer(droot[min(drank.id,na.rm = TRUE):nrow(droot),1])} #revised 20210626
    
    
    for(i in 1:length(drank))
    {
      if(length(bin.small.node)==0){ss=NA}else{ss=sapply(1:length(bs.path),function(j){if(drank[i] %in% bs.path[[j]]){which(bs.path[[j]]==drank[i])}else{NA}})}
      if(sum(!is.na(ss))>0)
      {
        if(length(bg.path)==0){sg=integer(0)}else{
          sg=sapply(1:length(bg.path),function(j){if(drank[i] %in% bg.path[[j]]){which(bg.path[[j]]==drank[i])}else{NA}})
        }
        if((sum(!is.na(ss))+sum(!is.na(sg)))>=2)
        {
          id.s=which(!is.na(ss))
          id.g=which(!is.na(sg))
          cl.s=sapply(1:length(id.s), function(u){bs.path.cum[[id.s[u]]][ss[id.s[u]]]})
          if(length(id.g)==0){cl.g=integer(0L)}else{cl.g=sapply(1:length(id.g), function(u){bg.path.cum[[id.g[u]]][sg[id.g[u]]]})}
          id.min=which.min(c(cl.g,cl.s))
          id.s.com=bin.small.node[id.s]
          if(id.min>length(cl.g))
          {
            sp.bin.temp[which(as.integer(sp.bin.temp) %in% as.integer(id.s.com))]=drank[i] #revised 20210626
            if(sum(sp.bin.temp==drank[i])<bin.size.limit)
            {
              if(length(id.g)==0)
              {
                bin.small.node=bin.small.node[-id.s]
                bs.path=bs.path[-id.s]
                bs.path.cum=bs.path.cum[-id.s]
                bs.path.len=bs.path.len[-id.s]
                bin.small.node[length(bin.small.node)+1]=drank[i]
                bs.path[[length(bs.path)+1]]=c(drank[i],node.path[[which(nodes==drank[i])]][[1]])
                bs.path.len[[length(bs.path.len)+1]]=c(0,node.path[[which(nodes==drank[i])]][[2]])
                bs.path.cum[[length(bs.path.cum)+1]]=cumsum(bs.path.len[[length(bs.path.cum)+1]])
              }else{
                id.min.g=id.g[which.min(cl.g)]
                sp.bin.temp[sp.bin.temp %in% drank[i]]=bin.good.node[id.min.g]
                bin.small.node=bin.small.node[-id.s]
                bs.path=bs.path[-id.s]
                bs.path.cum=bs.path.cum[-id.s]
                bs.path.len=bs.path.len[-id.s]
              }
            }else{
              bin.small.node=bin.small.node[-id.s]
              bs.path=bs.path[-id.s]
              bs.path.cum=bs.path.cum[-id.s]
              bs.path.len=bs.path.len[-id.s]
              bin.good.node=c(bin.good.node,drank[i])
              if(i==length(drank))
              {
                bg.path=c(bg.path,list(drank[i]));bg.path.len=c(bg.path.len,list(0));bg.path.cum=c(bg.path.cum,list(0))
              }else{
                bg.path=c(bg.path,list(c(drank[i],node.path[[which(nodes==drank[i])]][[1]])))
                bg.path.len=c(bg.path.len,list(c(0,node.path[[which(nodes==drank[i])]][[2]])))
                bg.path.cum=c(bg.path.cum,list(cumsum(bg.path.len[[length(bin.good.node)]])))
              }
            }
          }else{
            id.min.g=id.g[which.min(cl.g)]
            sp.bin.temp[sp.bin.temp %in% id.s.com]=bin.good.node[id.min.g]
            bin.small.node=bin.small.node[-id.s]
            bs.path=bs.path[-id.s]
            bs.path.cum=bs.path.cum[-id.s]
            bs.path.len=bs.path.len[-id.s]
          }
        }
      }
    }
    if(length(bin.small.node)>0){stop("something must be wrong. the small bins still exist. Stopped. ",date())}
    sp.bin[,2]=sp.bin.temp
    sp.bin[,3]=as.numeric(factor(sp.bin.temp))
    node.table.new=as.matrix(table(sp.bin.temp))
    sp.core.new=data.frame(bin.united.id.old=rownames(node.table.new),bin.united.tax.num=node.table.new,stringsAsFactors=FALSE)
    rownames(sp.core.new)=1:nrow(sp.core.new)
    newbin.sp.id=lapply(1:nrow(sp.core.new), function(i){rownames(sp.bin)[which(sp.bin[,3]==i)]})
    pd.state.new=sapply(1:nrow(sp.core.new),
                        function(i)
                        {
                          newpdid=match(newbin.sp.id[[i]],pd.spname)
                          pd.temp=stats::as.dist(pd[newpdid,newpdid])
                          c(max(pd.temp),mean(pd.temp),stats::sd(pd.temp))
                        })
    sp.core.new=cbind(sp.core.new,t(matrix(pd.state.new,nrow=3,ncol=nrow(sp.core.new))))
    colnames(sp.core.new)[3:5]=c("bin.pd.max","bin.pd.mean","bin.pd.sd")
  }else{
    sp.core.new=NA
    newbin.sp.id=NA
  }
  
  list(sp.bin=sp.bin,bin.united.sp=newbin.sp.id,bin.strict.sp=bin.sp.id,
       state.strict=sp.core,state.united=sp.core.new)
}