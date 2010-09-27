pcrProc <- function(filename)
{
  rawpcrdata <- data.frame(read.delim(filename, sep="\t", header=TRUE));
  pcrNorm <- sapply(colnames(rawpcrdata)[2:10], function(y) tapply(rawpcrdata[,y], rawpcrdata$X, function(x) mean(x)));
  pcrData <- t(pcrNorm);
  mycol <- c("red", "black");
  experiment <- c("wt", "coi1", "jin1", "jai3", "jut_", "jutOE", "DM");
  treatment <- c("MeJA", "mock");
  for(gene in rownames(pcrData))
  {
    par(mfrow=c(4,2));
    for(t in experiment)
    {
      x <- max(pcrData[gene,intersect(grep(treatment[1], colnames(pcrData)), grep(t, colnames(pcrData)))]) + 2;
      print(x);
      plot(pcrData[gene,intersect(grep(treatment[1], colnames(pcrData)), grep(t, colnames(pcrData)))], type="l", col="red", ylim=c(0,x), main=sprintf("%s - %st",gene, t), xlab="time steps", ylab="Relative expression");
      lines(pcrData[gene,intersect(grep(treatment[2], colnames(pcrData)), grep(t, colnames(pcrData)))]);
      legend("topright", legend=treatment, pch=15, pt.cex=1.0, col=mycol, ncol=2);
    }
    if (readline(sprintf("%s -- q to quit: ", gene)) == "q")
    {
      break;
    }
  }
 
 return(pcrData);
}


epsdevice <- function(epsFilename)
{
  postscript(epsFilename, width = 8, height = 6, paper = "special", onefile = FALSE, horizontal = FALSE);
}


createPlot <- function(filename, noise=0)
{
  plotName <- sprintf("jasmonate%d.eps", noise);
  d <- read.table(filename, sep="\t", header=T);
  epsdevice(plotName);
  par(cex.axis = 0.9);
  boxplot(fitness ~ model, data = d, notch = TRUE, xlab = " models", ylab = "optimised divergence", main = sprintf("Jasmonate models, noise = %d%s", noise, "%"), ylim = c(0, 35));
  dev.off()
}

getMatrixDistance <- function(data1, data2)
{
  if(length(setdiff(rownames(data1), rownames(data2))) != 0 ) 
  {
    print("error, datasets are not equal");
    break;
  }

  d <- 0;
  for(gene in rownames(data1))
  {
    if(sd(as.numeric(data1[gene,])) == 0 || sd(as.numeric(data2[gene,])) == 0  )
    {
      d = d + 1;
    }
    else 
    {
      d = d + (1 - cor(as.numeric(data1[gene,]), as.numeric(data2[gene,])));
    }
  }

  print("Matrix distance: ");
  print(d);

}
