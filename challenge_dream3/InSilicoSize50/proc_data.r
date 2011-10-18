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
  invisible(data);
}


readGoldStandardData <- function(filename)
{
  d <- read.table(filename, header = TRUE, sep = "\t");
  rownames(d) <- as.character(d[[1]]);
  d <- d[, 2:ncol(d)];
  return(d);
}

gsd <- readGoldStandardData("InSilicoSize50-Ecoli1-null-mutants.tsv");
processdata(gsd);
