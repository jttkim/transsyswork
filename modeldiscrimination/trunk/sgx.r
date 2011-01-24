readExpr <- function(exprName)
{
  d <- read.table(exprName, header = TRUE, sep = "\t");
  # make first column into rownames -- fix up routine that writes this sometime...
  if (colnames(d)[1] == "X")
  {
    rownames(d) <- as.character(d[, 1]);
    d <- d[, 2:ncol(d)];
  }
  return(d);
}


epsdevice <- function(epsFilename)
{
  postscript(epsFilename, width = 8, height = 6, paper = "special", onefile = FALSE, horizontal = FALSE);
}


plotProfile <- function(d, geneName)
{
  p <- as.numeric(d[geneName, ]);
  n <- c("-\n-", "-\nhk", "-\nc1", "-\nc2", "-\nc3", "+\n-", "+\nhk", "+\nc1", "+\nc2", "+\nc3");
  barplot(p, names = n, ylim = c(-8, 0), main = geneName);
}


profilePlots <- function(d)
{
  for (geneName in rownames(d))
  {
    epsdevice(sprintf("sgx_profile_%s.eps", geneName));
    plotProfile(d, geneName);
    dev.off();
  }
}
