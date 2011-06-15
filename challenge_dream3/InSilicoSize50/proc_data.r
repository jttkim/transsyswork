processdata <- function(data)
{
  for(i in 1:50)
  {
    rname <- paste("G", i, "(-/-)", sep="");
    t <- which(rownames(data) == rname);
    rownames(data)[t] <- paste("G", i, "KO", sep="");
  }

  data <- t(data);
  write.table(data, file="procdata.tsv", sep="\t", quote=FALSE);
}
