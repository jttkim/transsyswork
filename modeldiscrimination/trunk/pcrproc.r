pcrproc <- function(filename)
{
  rawpcrdata <- data.frame(read.delim(filename, sep="\t", header=TRUE));
  pcrNorm <- sapply(colnames(rawpcrdata)[2:10], function(y) tapply(rawpcrdata[,y], rawpcrdata$X, function(x) mean(x)));
  pcrData <- t(pcrNorm);
  mycol <- c("red", "black");
  experiment <- c("wt", "coi1", "jin1", "jai3", "jut_", "jutOE", "DM");
  treatment <- c("MJ", "mock");
  for(gene in rownames(pcrData))
  {
    par(mfrow=c(4,2));
    for(t in experiment)
    {
      x <- max(pcrData[gene,intersect(grep(treatment[1], colnames(pcrData)), grep(t, colnames(pcrData)))]) + 2;
      plot(pcrData[gene,intersect(grep(treatment[1], colnames(pcrData)), grep(t, colnames(pcrData)))], type="l", col="red", ylim=c(0,x), main=sprintf("%s - %st",gene, t), xlab="time steps", ylab="Relative expression");
      lines(pcrData[gene,intersect(grep(treatment[2], colnames(pcrData)), grep(t, colnames(pcrData)))]);
      legend("topright", legend=treatment, pch=15, pt.cex=1.0, col=mycol, ncol=2);
    }
    if (readline(sprintf("%s -- q to quit: ", gene)) == "q")
    {
      break;
    }
  }

}
