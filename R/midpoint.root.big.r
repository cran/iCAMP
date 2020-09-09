midpoint.root.big<-function (tree, pd.desc, pd.spname, pd.wd, nworker=4) 
{
  sp.name=tree$tip.label
  if(length(setdiff(sp.name,pd.spname))>0) stop("tree tip names and pd species names do not match.")
  
  requireNamespace("ape")
  requireNamespace("bigmemory")
  
  pd=bigmemory::attach.big.matrix(dget(paste0(pd.wd,"/",pd.desc)))
  maxdis=iCAMP::maxbigm(m.desc = pd.desc,m.wd = pd.wd,nworker = nworker,
                        rm.na = TRUE,size.limit = 10000*10000)
  
  dd <- maxdis$max.value
  ii <- as.vector(maxdis$row.col[1,])
  spp <- pd.spname[ii]
  nn <- which(tree$tip.label == spp[2])
  
  
  rerootn<-function (tree, node.number, position = NULL) 
  {
    # simplified from reroot in R package phytools 
    if (is.null(position)) 
      position <- tree$edge.length[which(tree$edge[, 2] == node.number)]
    
    splitTreen<-function (tree, split) 
    {
      # simplified from splitTree in R package phytools
      if (split$node > Ntip(tree))
      {
        tr2 <- ape::extract.clade(tree, node = split$node)
        tr2$root.edge <- tree$edge.length[which(tree$edge[, 2] == split$node)] - split$bp
        
        drop.claden<-function (tree, tip) 
        {
          nn <- if (!is.null(tree$node.label)){c(tree$node.label, "NA")}else{"NA"}
          tree <- ape::drop.tip(tree, tip, trim.internal = FALSE)
          while (sum(tree$tip.label %in% nn) > 1)
          {
            tree <- ape::drop.tip(tree,tree$tip.label[tree$tip.label %in% nn], trim.internal = FALSE)
          }
          tree
        }
        
        tr1 <- drop.claden(tree, tr2$tip.label)
        nn <- if (!is.null(tree$node.label)){c(tree$node.label, "NA")}else{"NA"}
        tr1$tip.label[which(tr1$tip.label %in% nn)] <- "NA"
        tr1$edge.length[match(which(tr1$tip.label == "NA"), tr1$edge[,2])] <- split$bp
      } else {
        tr2 <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = tree$tip.label[split$node], 
                    edge.length = tree$edge.length[which(tree$edge[,2] == split$node)] - split$bp, Nnode = 1L)
        class(tr2) <- "phylo"
        tr1 <- tree
        tr1$edge.length[match(which(tr1$tip.label == tr2$tip.label[1]),tr1$edge[, 2])] <- split$bp
        tr1$tip.label[which(tr1$tip.label == tr2$tip.label[1])] <- "NA"
      }
      trees <- list(tr1, tr2)
      class(trees) <- "multiPhylo"
      trees
    }
    tt <- splitTreen(tree, list(node = node.number, bp = position))
    p <- tt[[1]]
    d <- tt[[2]]
    tip <- if (length(which(p$tip.label == "NA")) > 0){"NA"}else{p$tip.label[which(p$tip.label %in% tree$node.label)]}
    p <- ape::root.phylo(p, outgroup = tip, resolve.root = TRUE)
    bb <- which(p$tip.label == tip)
    p$tip.label[bb] <- "NA"
    ee <- p$edge.length[which(p$edge[, 2] == bb)]
    p$edge.length[which(p$edge[, 2] == bb)] <- 0
    cc <- p$edge[which(p$edge[, 2] == bb), 1]
    dd <- setdiff(p$edge[which(p$edge[, 1] == cc), 2], bb)
    p$edge.length[which(p$edge[, 2] == dd)] <- p$edge.length[which(p$edge[,2] == dd)] + ee
    
    paste.treen<-function (tr1, tr2) 
    {
      # from paste.tree in R package phytool.
      if (length(tr2$tip) > 1)
      {
        temp <- tr2$root.edge
        tr2$root.edge <- NULL
        tr1$edge.length[match(which(tr1$tip.label == "NA"), tr1$edge[,2])] <- tr1$edge.length[match(which(tr1$tip.label =="NA"), tr1$edge[, 2])] + temp
      }
      tr.bound <- ape::bind.tree(tr1, tr2, where = which(tr1$tip.label == "NA"))
      return(tr.bound)
    }
    obj <- paste.treen(p, d)
    obj
  }
  tree <- rerootn(tree, nn, tree$edge.length[which(tree$edge[,2] == nn)])
  rootn=tree$edge[1,1]
  aa=which(tree$tip.label == spp[1])
  D<-0
  D.name<-aa
  while(aa!=rootn)
  {
    id=which(tree$edge[,2]==aa)
    Dn=tree$edge.length[id]+D[length(D)]
    D=c(D,Dn)
    aa=tree$edge[id,1]
    D.name=c(D.name,aa)
  }
  
  i <- 0
  while (D[i + 1] < (dd/2)) i <- i + 1
  tree <- rerootn(tree, as.numeric(D.name[i]), D[i + 1] - dd/2)
  list(tree=tree,max.pd=dd)
}