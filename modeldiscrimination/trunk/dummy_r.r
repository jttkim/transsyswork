plotdata  <- function(filename)
{

  data <- read.table(filename, header=TRUE, sep="\t")
  v=length(colnames(data))-2;
  par(mfrow=c(3,v))
  for(i in 1:length(levels(data$Array)))
    for (j in 2:(length(colnames(data))-1))
    {
      plot(data[which(data$Array == levels(data$Array)[i]),j], type='l', ylab=sprintf("%s",levels(data$Array)[i]), main=sprintf("%s",colnames(data)[j]), cex.lab=1.5, cex.main=1.5, col='red');
    }
}

plotjasmonate <- function(filename)
{

  label <- c("dummy","COI1","jasmonate","JIN1", "f_jin1", "JAZ"); 
  data <- read.table(filename, header=TRUE, sep="\t")
  v=length(colnames(label));
  par(mfrow=c(3,2))
  for(i in label)
  {
    plot(data[,i], type='l', ylab=sprintf("%s",i), main=sprintf("%s",i), cex.lab=1.5, cex.main=1.5, col='red');
  }
}


plotProfiles <- function(d)
{
  for (n in rownames(d))
  {
    p <- as.numeric(d[n, ]);
    names(p) <- colnames(d);
    barplot(p, main = n);
    if (readline(sprintf("%s -- q to quit: ", n)) == "q")
    {
      break;
    }
  }
}


traceColors <- function(nColors)
{
  n <- as.integer(ceiling(nColors / 3));
  s <- 1.0 / n;
  tc <- character();
  for (i in 0:(n - 1))
  {
    r <- as.integer(round((1.0 - i * s) * 255));
    g <- as.integer(round(i * s * 255));
    b <- as.integer(0);
    tc <- c(tc, sprintf("#%02x%02x%02x", r, g, b));
    r <- as.integer(0);
    g <- as.integer(round((1.0 - i * s) * 255));
    b <- as.integer(round(i * s * 255));
    tc <- c(tc, sprintf("#%02x%02x%02x", r, g, b));
    r <- as.integer(round(i * s * 255));
    g <- as.integer(0);
    b <- as.integer(round((1.0 - i * s) * 255));
    tc <- c(tc, sprintf("#%02x%02x%02x", r, g, b));
  }
  return(tc[1:nColors]);
}


plotTraces <- function(d, factorNameList = NULL)
{
  if (is.null(factorNameList))
  {
    factorNameList <- colnames(d)[2:ncol(d)];
  }
  arrayIndexList <- list();
  for (arrayName in as.character(unique(d[["Array"]])))
  {
    arrayIndexList[[arrayName]] <- which(d[["Array"]] == arrayName);
  }
  arrCol <- traceColors(length(arrayIndexList));
  names(arrCol) <- names(arrayIndexList);
  for (factorName in factorNameList)
  {
    factorTrace <- d[[factorName]];
    plot(1:length(factorTrace), factorTrace, type = "l", col = "white", main = factorName);
    for (arrayName in names(arrayIndexList))
    {
      lines(arrayIndexList[[arrayName]], factorTrace[arrayIndexList[[arrayName]]], type = "l", col = arrCol[arrayName]);
    }
    if (readline(sprintf("%s -- q to quit: ", factorName)) == "q")
    {
      break;
    }
  }
}


plotArrayTraces <- function(d, factorNameList = NULL, arrayNameList = NULL)
{
  if (is.null(factorNameList))
  {
    factorNameList <- colnames(d)[2:ncol(d)];
  }
  if (is.null(arrayNameList))
  {
    arrayNameList <- as.character(unique(d[["Array"]]));
  }
  for (factorName in factorNameList)
  {
    for (arrayName in arrayNameList)
    {
      plot(d[d[["Array"]] == arrayName, factorName], type = "l", main = sprintf("factor %s, array %s", factorName, arrayName));
      userQuit <- readline(sprintf("factor %s, array %s -- q to quit: ", factorName, arrayName)) == "q";
      if (userQuit)
      {
        break;
      }
    }
    if (userQuit)
    {
      break;
    }
  }
}


readSimsetExpression <- function(filename)
{
  d <- read.table(filename, header = TRUE, sep = "\t");
  rownames(d) <- as.character(d[[1]]);
  d <- d[, 2:ncol(d)];
  return(d);
}
