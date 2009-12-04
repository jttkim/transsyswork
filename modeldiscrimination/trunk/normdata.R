normdata<-function() {
 
}

getsequencemx <- function(filename, top, id) {
   
  values <- matrix(0,top,3)
  for(i in 1:top) {
    mx <-i*2
    if( i < 5 ) 
      name <- paste(filename,'0',mx,'01','_expr.txt',sep="")
    else
      name <- paste(filename,mx,'01','_expr.txt',sep="")
    values[i,] <- getvalues(name)
  }
 
  colnames(values) <- c("cv_all", "mean_FC", "sd_FC")
  va <- seq(2,(top*2), by=2)
  values <- data.frame(values,row.names=va)
  filename <- paste("mxnoise1id", id,".txt",sep="")
  write.table(values, file=filename, sep='\t')
  return(values)

}


getsequencegn <- function(filename, min,max, id) {
   
  values <- matrix(0,((max-min)+1),3)
  for(i in min:max) {
    if( i < 5 ) 
      name <- paste(filename,'0',mx,'01','_expr.txt',sep="")
    else
      name <- paste(filename,'n',i,'e',(i*id),'0101','_expr.txt',sep="")
    values[(i-14),] <- getvalues(name)
  }
 
  colnames(values) <- c("cv_all", "mean_FC", "sd_FC")
  va <- seq(min,max, by=1)
  values <- data.frame(values,row.names=va)
  print(values)
  filename <- paste("tpnoise1", id,".txt",sep="")
  write.table(values, file=filename, sep='\t')
  return(values)

}


getsequencetop <- function(filename, n, id, t, p, net) {
  
  inx <- 1   
  values <- matrix(0,(t*p),5)
  for(i in 1:t) {
    for(j in 1:p) {
      if( i < 10 ) 
        name <- paste(filename,'n',n,'e',(n*id),'0',i,'0',j,'_expr.txt',sep="")
      else
        name <- paste(filename,'n',n,'e',(n*id),i,'0',j,'_expr.txt',sep="")
      values[inx,] <- c(i,j,getvalues(name))
      inx <- inx + 1
    }
  }
 
  colnames(values) <- c("topology", "parameter", "cv_all", "lr_mean", "lr_sd")
  return(values)
  

}

getvalues <- function(filename){

  data<-read.table(filename, sep='\t',header=TRUE, row.names=1)
  cv<-sd(as.vector(as.matrix(data))) / mean(as.vector(as.matrix(data)))
  
  data<-shiftdata(data)
  data<-normalise2mean(data)

  for(i in 2:ncol(data)) 
      data[,i]<-log2(data[,i]/data[,1])

  mean<- mean(as.vector(as.matrix(data[,2:ncol(data)])))
  sd<-sd(as.vector(as.matrix(data[2:ncol(data)])))

  return(c(cv, mean, sd))
}

normalise2mean <- function(data) {

  for(i in 1:ncol(data)) 
    data[,i]<-data[,i]/mean(as.vector(as.matrix(data[,i]))) # Per array normalisation
  return(data)

}

shiftdata <- function(data) {

  d<-as.vector(as.matrix(data[,1:ncol(data)]))
  d<-as.vector(d)
  logratio_offset<- sd(d) * 0.01
  m <- min(data)
  mfloor<-logratio_offset-m
  data<-data+mfloor 
  
  return (data)

}

plotdata <- function(data) {

  par(mfrow=c(5,2))
  for (i in 1:nrow(data)) {
    v<-as.list(as.vector(strsplit(rownames(data)[i], NULL)))
    index<-(which((colnames(data) == paste("a",v[[1]][4],v[[1]][5],sep="")))*-1)
    plot(as.vector(as.matrix(data[i,c(index)])), col="red", axes=FALSE, xlab="Gene mutant", ylab="Av. Expres. Level")
    #plot(as.vector(as.matrix(data[i,])), col="red", axes=FALSE, xlab="Gene mutant", ylab="Av. Expres. Level")
    axis(1, 1:ncol(data), colnames(data))
    axis(2)
    box()
  }

}

plotbox <- function(mydata, n, id, noise, net, sepa) {
  
  if (sepa == 't')
    data <- read.table(mydata, header=TRUE, sep='\t')
  else
    data <- read.table(mydata, header=TRUE, sep='')
  colnames(data) <- c("rw_operation", "rw_repetition","fitness")
  m <- n*id
  #xxx <-bquote("Net "~.(net)~ " n="~.(n)~ " m="~.(m))
  xxx <-bquote("Net "~.(net)~ " n="~.(n)~ " m="~.(m)~"  noise="~.(noise)~"%")
  #boxplot(fitness ~ rw_operation, data=data, main=xxx, col="green", notch=TRUE, xlab="Rewiring operations", ylab="Repetitions", sub="1% noise added, global average, logratio, amx,rmx=10")
  boxplot(fitness ~ rw_operation, data=data, main=xxx, col="grey", notch=TRUE, xlab="Rewiring operations", ylab="score", sub="(b)")

}


plotplot <-function( er, pl) {
  #pdf("er_pl.pdf", width=8.5, height=11, paper="a4")
  #jpeg("sd.jpeg", width=8.5, height=11, paper="a4")
  #par(mfrow=c(2,2))
  #boxplot(lr_sd ~ topology, data = er, col = "red", main="Net ER, n=30, m=60", xlab="Topology", ylab="lr sd", notches=T)
  #boxplot(lr_sd ~ topology, data = pl, col = "blue", main="Net Powerlaw, n=30, m=60", xlab="Topology", ylab="lr sd", notches=T)
  plot(as.vector(er[,5]), col="red", main="lr sd, ER vs Powerlaw", ylab="lr sd", ylim=c(0,1.5))
  points(as.vector(pl[,5]), col="blue" )
  legend("topright", legend=c("ER","PL"), cex=1.3, bty="n", pch=1, pt.cex=1.8, col=c("red", "blue"))
  abline(h=0.44, col="red")
  #dev.off()
}


plot_utest <-function(data) {

  data <-read.delim(data, sep='', header=T)
  mtag <-c(1, 2, 3, 4, 5, 6, 7, 9, 11, 13, 15, 18, 22, 27, 32, 38, 46, 55, 66)
  mpvalue <- matrix(0, (length(mtag)+171),2)
  colnames(mpvalue)<- c("rewiring", "p-value")
  for(i in 1:length(mtag)) {
    result<-wilcox.test(data[which(data$rewiring == 0),3], data[which(data$rewiring == mtag[i]),3], paired=FALSE, conf.level = 0.95)
    mpvalue[i,1:2]<-c(mtag[i],result[[3]]*10)
  }

  ml <- length(mtag)
  for(i in 1:(length(mtag)-1)) {
    for(j in (i+1):length(mtag)) {
      result<-wilcox.test(data[which(data$rewiring == mtag[i]),3], data[which(data$rewiring == mtag[j]),3], paired=FALSE, conf.level = 0.95)
      comb<-paste(mtag[i],"-",mtag[j],sep='')
      ml <- ml + 1
      mpvalue[ml,1:2]<-c(comb,result[[3]])
    }
  }

  #plot(mpvalue[2:length(mtag),2], pch=0.5, axes=FALSE, xlab='rewiring operations', ylab='p-value')
  #axis(1, 1:length(mtag), mtag[1:length(mtag)])
  #axis(2)
  #box()
  return(mpvalue)
}
