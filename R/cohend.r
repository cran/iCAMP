cohend<-function(treat,control,paired=FALSE)
{
  if(paired)
  {
    if(length(treat)!=length(control)){stop("cannot do paired test for unequal sample size.")}
    x1=as.numeric(treat)-as.numeric(control)
    x2=0
    sd1=stats::sd(x1)
    cod=mean(x1)/sd1
    n1=length(x1)
    n2=1
  }else if(length(treat)==1|length(control)==1){
    x1=as.numeric(treat)
    x2=as.numeric(control)
    n1=length(x1)
    n2=length(x2)
    if(n1==1){sd1=0}else{sd1=stats::sd(x1)}
    if(n2==1){sd2=0}else{sd2=stats::sd(x2)}
    sd12=((((n1-1)*(sd1^2))+((n2-1)*(sd2^2)))/(n1+n2-2))^0.5
    cod=(mean(x1)-mean(x2))/sd12
  }else{
    x1=as.numeric(treat)
    x2=as.numeric(control)
    n1=length(x1)
    n2=length(x2)
    sd12=((((n1-1)*(stats::sd(x1)^2))+((n2-1)*(stats::sd(x2)^2)))/(n1+n2-2))^0.5
    cod=(mean(x1)-mean(x2))/sd12
  }
  if(is.na(cod))
  {
    output=list(d=0,sd=NA,magnitude="negligible",paired=paired)
  }else{
    var.cod=(((n1+n2)/(n1*n2))+((cod^2)/(2*(n1+n2-2))))*((n1+n2)/(n1+n2-2))
    sd.cod=var.cod^0.5
    
    abs.d=abs(cod)
    if(abs.d<0.2){mag="negligible"}else if(abs.d<0.5){mag="small"}else if(abs.d<0.8){mag="medium"}else {mag="large"}
    output=list(d=cod,sd=sd.cod,magnitude=mag,paired=paired)
  }
  output
}