change.sigindex<-function(icamp.output,sig.index=c("Confidence","SES.RC","SES","RC"),
                          detail.save=TRUE,detail.null=FALSE,
                          ses.cut = 1.96,rc.cut = 0.95,conf.cut=0.975)
{
  sig.index=sig.index[1]
  if(!(sig.index %in% c("SES.RC","Confidence","RC","SES"))){stop("wrong sig.index for change.sigindex")}
  
  sig.index.old=icamp.output$detail$setting$sig.index
  ses.cut.old=icamp.output$detail$setting$ses.cut
  rc.cut.old=icamp.output$detail$setting$rc.cut
  conf.cut.old=icamp.output$detail$setting$conf.cut
  
  if(sig.index.old==sig.index & ses.cut.old==ses.cut & rc.cut.old==rc.cut & conf.cut.old==conf.cut)
  {
    warning("Nothing need to change, return original results.")
    output=icamp.output
  }else{
    if(sig.index=="SES.RC")
    {
      sig.ind.phy="SES";sig.ind.tax="RC"
    }else if(sig.index=="SES"){
      sig.ind.phy<-sig.ind.tax<-"SES"
    }else if(sig.index=="RC"){
      sig.ind.phy<-sig.ind.tax<-"RC"
    }else if(sig.index=="Confidence"){
      sig.ind.phy<-sig.ind.tax<-"Confidence"
    }else{stop("wrong sig.index in change.sigindex.")}
    
    detail=icamp.output$detail
    detail$setting$sig.index=sig.index
    detail$setting$detail.save=detail.save
    detail$setting$detail.null=detail.null
    detail$setting$ses.cut=ses.cut
    detail$setting$rc.cut=rc.cut
    detail$setting$conf.cut=conf.cut
    
    sigindexn<-function(obs,rands,special3=NULL,sig.ind,obj.name)
    {
      randar=sapply(1:length(rands),function(i){as.matrix(rands[[i]][,3:ncol(rands[[i]]),drop=FALSE])},simplify = "array")
      obsa=sapply(1:length(obs),function(i){iCAMP::dist.3col(obs[[i]])[,3]})
      
      if(sig.ind=="RC")
      {
        obsar=sapply(1:ncol(obsa),function(i){matrix(obsa[,i],nrow=nrow(randar),ncol=ncol(randar))},simplify = "array")
        alpha1=apply(obsar==randar,c(1,3),sum)
        alpha=apply(obsar>randar,c(1,3),sum)
        alpha=(alpha+0.5*alpha1)/(dim(randar)[2])
        out1=2*alpha-1
        if(!is.null(special3)){special=special3$special.rc}
      }else if(sig.ind=="SES"){
        out1=(obsa-apply(randar,c(1,3),mean))/(apply(randar,c(1,3),stats::sd))
        if(!is.null(special3)){special=special3$special.ses}
      }else if(sig.ind=="Confidence"){
        obsar=sapply(1:ncol(obsa),function(i){matrix(obsa[,i],nrow=nrow(randar),ncol=ncol(randar))},simplify = "array")
        alpha=(apply(obsar>randar,c(1,3),sum))/(dim(randar)[2])
        alpha2=(apply(obsar<randar,c(1,3),sum))/(dim(randar)[2])
        alpha[which(alpha2>alpha, arr.ind = TRUE)]=-alpha2[which(alpha2>alpha, arr.ind = TRUE)]
        out1=alpha
        if(!is.null(special3)){special=special3$special.conf}
      }
      colnames(out1)=paste0(obj.name,".bin",1:ncol(out1))
      if(!is.null(special3))
      {
        specialm=special[,3:ncol(special)]
        if(sum(specialm!=0)>0)
        {
          out1[which(specialm!=0,arr.ind = TRUE)]=specialm[which(specialm!=0,arr.ind = TRUE)]
        }
      }
      out1[is.na(out1)]=0
      data.frame(rands[[1]][,1:2],out1,stringsAsFactors = FALSE)
    }
    
    objnameck<-function(indname,sig.ind.phy,sig.ind.tax)
    {
      if(grepl("BC",indname) | grepl("bray",tolower(indname)))
      {
        if(sig.ind.tax=="Confidence"){outk=paste0("CBray",substring(indname,nchar(indname)))}else{
          outk=paste0(sig.ind.tax,"bray",substring(indname,nchar(indname)))
        }
      }else{
        if(sig.ind.phy=="SES")
        {
          if(grepl("bMPD",indname)){outk=sub("bMPD","bNRI",indname)}
          if(grepl("bMNTD",indname)){outk=sub("bMNTD","bNTI",indname)}
        }else if(sig.ind.phy=="Confidence"){
          outk=paste0("C",indname)
        }else{
          outk=paste0("RC",indname)
        }
      }
      outk
    }
    
    idrand=grep(".rand$",names(icamp.output$rand))
    for(i in 1:length(idrand))
    {
      indxi=sub(".rand","",names(icamp.output$rand)[idrand[i]])
      message("Now re-calculate i=",i," in ",length(idrand),", ",indxi,". ",date())
      randxi=icamp.output$rand[[idrand[i]]]
      idobsi=grep(paste0(substr(indxi,1,nchar(indxi)-1),".obs"),names(icamp.output$rand))[1]
      obsxi=icamp.output$rand[[idobsi]]
      idspeci=grep(paste0("Sig",indxi),names(icamp.output$special.crct))[1]
      specialxi=icamp.output$special.crct[[idspeci]]
      objnamexi=objnameck(indname = indxi, sig.ind.phy = sig.ind.phy, sig.ind.tax = sig.ind.tax)
      iddet=grep(paste0("Sig",indxi),names(detail))[1]
      detail[[iddet]]=sigindexn(obs=obsxi,rands=randxi,special3 = specialxi,sig.ind=sig.ind.phy,obj.name=objnamexi)
    }
    
    if(sig.ind.phy=="SES"){sig.phy.cut<-sig.phy2.cut<-ses.cut}else if(sig.ind.phy=="RC"){sig.phy.cut<-sig.phy2.cut<-rc.cut}else{sig.phy.cut<-sig.phy2.cut<-conf.cut}
    if(sig.ind.tax=="SES"){sig.tax.cut<-ses.cut}else if(sig.ind.tax=="RC"){sig.tax.cut<-rc.cut}else{sig.tax.cut<-conf.cut}  
    bin.weight=icamp.output$detail$bin.weight
    
    phy1.ids=grep("bMPD",names(detail))
    phy2.ids=grep("bMNTD",names(detail))
    phy.ids=c(phy1.ids,phy2.ids)
    tax.ids=grep("BC",names(detail))
    
    qpres=list()
    for(i in 1:length(phy.ids))
    {
      phynamei=sub(".bin1","",colnames(detail[[phy.ids[i]]])[3])
      for(j in 1:length(tax.ids))
      {
        taxnamej=sub(".bin1","",colnames(detail[[tax.ids[j]]])[3])
        message("Now quantify process importance based on ",phynamei," and ",taxnamej,". ",date())
        qpij=iCAMP::qp.bin.js(sig.phy.bin=detail[[phy.ids[i]]], sig.phy2.bin=NULL, 
                              sig.tax.bin=detail[[tax.ids[j]]],
                              bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                              sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                              check.name=FALSE)
        qpres=c(qpres,list(qpij))
        names(qpres)[length(qpres)]=paste0(phynamei,taxnamej)
      }
    }
    
    if((length(phy1.ids)*length(phy2.ids))>0)
    {
      for(i in 1:length(phy1.ids))
      {
        phynamei=sub(".bin1","",colnames(detail[[phy1.ids[i]]])[3])
        iora=substring(phynamei,nchar(phynamei))
        jids=grep(paste0(iora,"$"),names(detail)[phy2.ids])
        if(length(jids)>0)
        {
          for(j in jids)
          {
            phynamej=sub(".bin1","",colnames(detail[[phy2.ids[j]]])[3])
            for(k in 1:length(tax.ids))
            {
              taxnamek=sub(".bin1","",colnames(detail[[tax.ids[k]]])[3])
              message("Now quantify process importance based on ",phynamei,", ",phynamej,", and ",taxnamek,". ",date())
              qpijk=iCAMP::qp.bin.js(sig.phy.bin=detail[[phy1.ids[i]]],
                                     sig.phy2.bin=detail[[phy2.ids[j]]], 
                                     sig.tax.bin=detail[[tax.ids[k]]],
                                     bin.weight=bin.weight, sig.phy.cut=sig.phy.cut, 
                                     sig.phy2.cut=sig.phy2.cut, sig.tax.cut=sig.tax.cut,
                                     check.name=FALSE)
              qpres=c(qpres,list(qpijk))
              names(qpres)[length(qpres)]=paste0(phynamei,phynamej,taxnamek)
            }
          }
        }
      }
    }
    
    output=qpres
    if(detail.save)
    {
      detail$processes=qpres
      output=c(output,list(detail=detail))
    }
    if(detail.null)
    {
      output=c(output,list(rand=icamp.output$rand,special.crct=icamp.output$special.crct))
    }
  }
  output
}