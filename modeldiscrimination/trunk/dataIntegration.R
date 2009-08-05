library(affy)
library(ArrayExpress)


dataIntegration  <- function() { 


}


dataint <- function() {

  array_names <- c("E-GEOD-10646","E-GEOD-11216","E-GEOD-2169","E-GEOD-6151","E-GEOD-9702","E-GEOD-9955","E-MEXP-1787","E-MEXP-883")
  #array_names <- c("E-GEOD-10646","E-GEOD-11216")
  system("clear")
  cv_array <- c() 

  for(i in 1:length(array_names)) { 
    label <- c() 
    eset <- ArrayExpress(array_names[i])
    name <- rma(eset)
    mydata <- exprs(name)
    colnames(mydata) <- phenoData(name)$Hybridization.Name
    mydata <- deleteaffx(mydata)
    mydata <- 2^mydata

    if(i==1) {
      meta <- data.frame(mydata[,1])
      rownames(meta) <- rownames(mydata)
      colnames(meta) <- c("temp")
      lr_meta <- meta
      fc_meta <- meta
      nk_meta <- meta
    }

    lp <- dim(mydata)[2]
    system("clear")
    print(phenoData(name)$Hybridization.Name)
  
    an=0
    while (an != 1) {
      print("Entry knockout gene names. Press enter to finish" )
      nkname  <-  scan(what=character(),sep=",",quiet=T,nlines=1)
      print("Is this correct? :")
      print(nkname)
      print("1: correct, 0: otherwise")
      an <- scan(what=integer(),nmax=1)
    }
  
    an=0
    while (an != 1) {
      print("What column numbers are replicates of wt? e.g.(1,2,3). Press enter to finish")
      wt <- scan(what=numeric(), sep="," ,quiet=T, nlines=1)
      print("Is this correct? :")
      print(phenoData(name)$Hybridization.Name[wt])
      print("1: correct, 0: otherwise")
      an <- scan(what=integer(),nmax=1)
    }
    mydata <- data.frame(mydata, label=apply(mydata[,wt], 1, mean))

    nk_meta <- data.frame(merge(nk_meta,mydata[,wt],by.x=0, by.y=0),row.names=1)
    cv_array <- append(cv_array,mean(apply(mydata[,wt],1,"sd"))/mean(mydata[,dim(mydata)[2]]))
    label <- append(label,(paste('wt','_avg',sep='')))
    for(j in 1:length(nkname)) {
      an=0
      while (an != 1) {
        al <- searchcolumn(phenoData(name)$Hybridization.Name,nkname[j])
        ctitle <- bquote("Which colnumbers are (mock) replicates of "~.(nkname[j])~" Press enter to finish")
        print(ctitle)
	print(al)
        nk <- scan(what=numeric(),sep=",",quiet=T,nlines=1)
        print("Is this correct? :")
        print(phenoData(name)$Hybridization.Name[nk])
        print("1: correct, 0: otherwise")
        an <- scan(what=integer(),nmax=1)
      }
      label <- append(label,(paste(nkname[j],'_avg',sep='')))
      mydata <- data.frame(mydata, label=apply(mydata[,nk], 1, mean))
      nk_meta <- data.frame(merge(nk_meta,mydata[,nk],by.x=0, by.y=0),row.names=1)
      cv_array <- append(cv_array,mean(apply(mydata[,nk],1,"sd"))/mean(mydata[,dim(mydata)[2]]))
    }

    print(mydata[1,])
    mydata <- mydata[,c((lp+1):dim(mydata)[2])]
    colnames(mydata) <- label
    meta <- data.frame(merge(meta,mydata,by.x=0, by.y=0),row.names=1)
    
    mydata <- normalise2mean(mydata)

    lp <- dim(mydata)[2]
    for(k in 1:length(nkname)) {
      mydata <- data.frame(mydata, ratio=log2(mydata[,k+1]/mydata[,1]))
    }
    if(length(nkname)==1) {
    
      mydata <- mydata[,c(lp:dim(mydata)[2])]
      colnames(mydata)[1] <- c("temp")
      colnames(mydata)[2] <- nkname[1]
    }
    else {
      mydata <- mydata[,c((lp+1):dim(mydata)[2])]
      colnames(mydata) <- nkname
    }
    print(mydata[1,]) 
    lr_meta <- data.frame(merge(lr_meta,mydata,by.x=0, by.y=0),row.names=1)
    system("clear")
    rm(mydata,name,lp,nkname,label,j,k)
  }

  a <- deletecolumn(colnames(lr_meta))
  a <- a*(-1)
  lr_meta <- lr_meta[,a]
  write.table(lr_meta,file="lrmeta.txt",sep='\t') 

  a <- deletecolumn(colnames(meta))
  a <- a*(-1)
  meta <- meta[,a]
  write.table(meta,file="avmeta.txt",sep='\t') 

  a <- deletecolumn(colnames(nk_meta))
  a <- a*(-1)
  nk_meta <- nk_meta[,a]
  write.table(nk_meta,file="nkmeta.txt",sep='\t') 

  # Calculate global coefficient of variability (across expression matrix)
  cv <- as.numeric(c(statistics(meta)))
  print(cv[2]/cv[1])

  return(cv_array)
}


statistics <- function(data) {

  c <- dim(data)[2]
  sd <- sd(as.vector(as.matrix(data)[,1:c]))
  av <- mean(as.vector(as.matrix(data)[,1:c]))
  
  return(av,sd)

}


deletecolumn <- function(array_col) {

  a <- c()
  for (i in 1:length(array_col)) { 
    if (charmatch("temp", array_col[i], nomatch=0) != 0) 
      a <- append(a,i)
  }
  return(a)
}


searchcolumn <- function(array_label,gene_name) {
  
  a <- c()
  for(i in 1:length(array_label)) {
    y <- array_label[i]
    if(length(grep(gene_name, as.character(y))) !=0 ) {
      name <- paste(y," colnumber: ",i,sep="")
      a <- append(a,name)
    }
  }
  return(a)
}


getfc <- function(mydata) {

  mydata <- data.frame(mydata, FC_mkk1=log2(mydata[,2]/mydata[,1]))
  mydata <- data.frame(mydata, FC_mkk2=log2(mydata[,3]/mydata[,1]))
  mydata <- data.frame(mydata, FC_arf2=log2(mydata[,5]/mydata[,4]))
  mydata <- data.frame(mydata, FC_rre1=log2(mydata[,7]/mydata[,6]))
  mydata <- data.frame(mydata, FC_rre2=log2(mydata[,8]/mydata[,6]))
  mydata <- data.frame(mydata, FC_aih=log2(mydata[,10]/mydata[,9]))
  mydata <- data.frame(mydata, FC_aah=log2(mydata[,11]/mydata[,9]))
  mydata <- data.frame(mydata, FC_arh=log2(mydata[,12]/mydata[,9]))
  mydata <- data.frame(mydata, FC_wlh=log2(mydata[,13]/mydata[,9]))
  mydata <- data.frame(mydata, FC_mxh=log2(mydata[,14]/mydata[,9]))
  mydata <- data.frame(mydata, FC_dex=log2(mydata[,16]/mydata[,15]))
  mydata <- data.frame(mydata, FC_mil4=log2(mydata[,18]/mydata[,17]))
  mydata <- data.frame(mydata, FC_sid2=log2(mydata[,19]/mydata[,17]))
  mydata <- data.frame(mydata, FC_jin1=log2(mydata[,21]/mydata[,20]))
  mydata <- data.frame(mydata, FC_ccr1=log2(mydata[,23]/mydata[,22]))

}


normalise2mean <- function(data) {

  for(i in 1:ncol(data)) 
    data[,i] <- data[,i]/mean(as.vector(as.matrix(data[,i])))

  return(data)

}

shiftdata <- function(data) {

  d <- as.vector(as.matrix(data[,1:ncol(data)]))
  d <- as.vector(d)
  logratio_offset <- sd(d) * 0.01
  m <- min(data)
  mfloor <- logratio_offset-m
  data <- data+mfloor

  return(data)

}


deleteaffx <- function(data) {

  a <- c()
  for (i in 1:nrow(data)){
    if (charmatch("AFFX", rownames(data)[i], nomatch=0) == 0)
      a <- append(a,i)
  }
  data <- data.frame(data[a,])

  return(data)

}


affy2ath <- function(name) {

  library("ath1121501.db")
  affy_sample <- name
  test_sample <- as.vector(unlist(mget(affy_sample, ath1121501ACCNUM, ifnotfound=NA)))

  return(test_sample)

}

