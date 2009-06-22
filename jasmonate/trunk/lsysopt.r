library(jtkstuff);
library(transsys);


getSymbolRecords <- function(trlName, transsysName, numTimesteps)
{
  cmd <- sprintf("ltransexpr -t %s -n %d %s", transsysName, as.integer(numTimesteps), trlName);
  print(cmd);
  p <- pipe(cmd);
  d <- read.table(p, header = TRUE);
  # close(p);
  return(d);
}


getDisparitySets <- function(r, ruleName, factorName)
{
  a <- r[["ruleName"]] == ruleName;
  s <- list();
  s[["active"]] <- r[[factorName]][a];
  s[["inactive"]] <- r[[factorName]][!a];
  attr(s, "ruleName") <- ruleName;
  attr(s, "factorName") <- factorName;
  return(s);
}


stddevDisparity <- function(r, ruleName, factorName)
{
  s <- getDisparitySets(r, ruleName, factorName);
  ma <- mean(s[["active"]]);
  sa <- sd(s[["active"]]);
  mi <- mean(s[["inactive"]]);
  si <- sd(s[["inactive"]]);
  print(sprintf("rule %s, factor %s: act = %e +- %e, inact = %e +- %e", ruleName, factorName, ma, sa, mi, si));
  mdiff <- abs(ma) - abs(mi);
  if (mdiff == 0.0)
  {
    sddisp <- 1.0;
  }
  else
  {
    sm <- (sa + si) / abs(ma - mi);
    sddisp <- sm / (1.0 + sm);
  }
  return(sddisp);
}


getFinalStats <- function(d, colNames = NULL)
{
  if (is.null(colNames))
  {
    colNames <- colnames(d);
    if ("tp_index" %in% colNames)
    {
      colNames = colNames[colNames != "tp_index"];
    }
    if ("restart_index" %in% colNames)
    {
      colNames = colNames[colNames != "restart_index"];
    }
  }
  dfinal <- data.frame("tp_index" = numeric(), "restart_index" = numeric());
  for (n in colNames)
  {
    dfinal[[n]] <- numeric();
  }
  for (tpIndex in unique(d[["tp_index"]]))
  {
    tpLogical <- d[["tp_index"]] == tpIndex;
    for (restartIndex in unique(d[["restart_index"]]))
    {
      restartLogical <- d[["restart_index"]] == restartIndex;
      if (any(tpLogical & restartLogical))
      {
        i <- max(which(tpLogical & restartLogical));
        ## print(sprintf("tp_index = %d, restart_index = %d, i = %d", as.integer(tpIndex), as.integer(restartIndex), as.integer(i)));
        dfinal[nrow(dfinal) + 1, ] <- c(tpIndex, restartIndex, d[i, colNames]);
      }
    }
  }
  return(dfinal);
}


plotCurves <- function(d, splitFactor, colNames = NULL, ...)
{
  splitFactor <- as.factor(splitFactor);
  if (is.null(colNames))
  {
    colNames <- colnames(d);
  }
  for (n in colNames)
  {
    plotListRainbow(tapply(d[[n]], splitFactor, function(x) {x;}), main = n, ylab = n, ...);
    invisible(readline(sprintf("%s -- hit return", n)));
  }
}
