source("jasmonate.r");


epsdevice <- function(epsname)
{
  postscript(epsname, width=8, height=6, paper="special", onefile=FALSE, horizontal=FALSE);
  par(cex = 1.5, cex.main = 1.3);
}


integratedPanelsEps <- function(data, geneGroup, epsNameFormat)
{
  for (groupName in names(geneGroup))
  {
    epsFileName <- sprintf(epsNameFormat, groupName);
    epsdevice(epsFileName);
    par(cex.axis = 0.5);
    # print(lineProfiles);
    plotHighlightPanels(integrateControls(data), geneGroup[groupName], lineProfiles, xstyle = "stacked", waitFunc = print);
    dev.off();
  }
}


profilePlotEps <- function(baseName, geneNames, ...)
{
  profileFile <- sprintf("%s.dat", baseName);
  pp <- read.table(profileFile, header = TRUE);
  for (geneName in geneNames)
  {
    epsFileName <- sprintf("%s_%s.eps", baseName, geneName);
    epsdevice(epsFileName);
    plotProfileSets(pp, geneName, waitFunc = function(x) { print(x); });
    dev.off();
  }
}


optimisationPlotEps <- function(baseFormat)
{
  fnameFormat <- sprintf("%s_log.dat.gz", baseFormat);
  print(fnameFormat);
  l <- readOptimisationLogs(fnameFormat, fileOpenFunc = gzfile);
  epsFileName <- sprintf("%s_optbox.eps", sprintf(baseFormat, ""));
  epsdevice(epsFileName);
  par(cex.axis = 0.5);
  optimisationBoxplot(l, notch = TRUE);
  dev.off();
}


optimisationPlotEpsGrad <- function(baseFormat, nmut)
{
  fnameFormat <- sprintf("%s_log.dat.gz", baseFormat);
  print(fnameFormat);
  # l <- readOptimisationLogs(fnameFormat, fileOpenFunc = gzfile);
  l <- readFrames(sprintf("%02d", nmut), fnameFormat, gzfile)
  epsFileName <- sprintf("%s_optbox.eps", sprintf(baseFormat, ""));
  print(epsFileName);
  epsdevice(epsFileName);
  par(cex.axis = 0.5);
  optimisationBoxplot(l, notch = TRUE);
  dev.off();
}


optimisationPlotEpsGradPolished <- function(objectiveSpec, nmut, rndseedList, ylab)
{
  finalObjList <- readAllGradFinalStats(objectiveSpec, nmut, rndseedList);
  epsFileName <- sprintf("newphytologist_testjfullmut_grad_tfhypplus_%s_optbox.eps", objectiveSpec);
  print(epsFileName);
  epsdevice(epsFileName);
  par(cex.axis = 0.9, cex.lab = 1.4, mar = c(5, 5, 1, 1) + 0.1);
  boxplot(finalObjList, notch = TRUE, xlab = "number of perturbations", ylab = ylab);
  dev.off();
}


readAllGradFinalStats <- function(objectiveSpec, nmut, rndseedList)
{
  finalObjList <- NULL;
  rndseedCol <- integer(0);
  for (rndseed in rndseedList)
  {
    # testjfullmut33_010_grad_tfhypplus_squaresum_logratio_log.dat.gz
    fnameFormat <- sprintf("testjfullmut%%s_%03d_grad_tfhypplus_%s_log.dat.gz", as.integer(rndseed), objectiveSpec);
    l <- readFrames(sprintf("%02d", nmut), fnameFormat, gzfile);
    f <- lapply(l, function(x) {getFinalStats(x)[["obj"]];});
    numRestarts <- max(unlist(lapply(f, length)));
    # print(numRestarts);
    if (numRestarts != min(unlist(lapply(f, length))))
    {
      stop("variation in number of restarts");
    }
    rndseedCol <- c(rndseedCol, as.integer(rep(rndseed, numRestarts)));
    if (is.null(finalObjList))
    {
      finalObjList <- f;
    }
    else
    {
      for (n in names(finalObjList))
      {
        finalObjList[[n]] <- c(finalObjList[[n]], f[[n]]);
      }
    }
  }
  attr(finalObjList, "rndseedColumn") <- rndseedCol;
  return(finalObjList);
}


optimisationPlotGradAll <- function(objectiveSpec, nmut, rndseedList)
{
  finalObjList <- readAllGradFinalStats(objectiveSpec, nmut, rndseedList);
  x <- numeric(0);
  y <- numeric(0);
  for (n in names(finalObjList))
  {
    x <- c(x, rep(as.integer(n), length(finalObjList[[n]])));
    y <- c(y, finalObjList[[n]]);
  }
  plot(x, y);
  readline("Hit return");
  rc <- rainbow(length(rndseedList));
  for (i in seq(along = rndseedList))
  {
    points(x[attr(finalObjList, "rndseedColumn") == rndseedList[i]], y[attr(finalObjList, "rndseedColumn") == rndseedList[i]], col = rc[i]);
  }
  readline("Hit return");
  boxplot(finalObjList, notch = TRUE);
}


nmut <- c(0, 2, 4, 6, 8, 10, 12, 14, 17, 20, 24, 28, 33);
rndseedList <- 1:10;
optimisationPlotEpsGradPolished("correlation", nmut, rndseedList, expression(D[corr]));
optimisationPlotEpsGradPolished("correlation_logratio", nmut, rndseedList, expression(D[corr]^"(l)"));
optimisationPlotEpsGradPolished("squaresum", nmut, rndseedList, expression(D[ssq]));
optimisationPlotEpsGradPolished("squaresum_logratio", nmut, rndseedList, expression(D[ssq]^"(l)"));

# x <- readFrames(sprintf("%02d", nmut), "testjfullmut%s_001_grad_correlation_tfhpplus_r5_log.dat.gz", gzfile)

# optimisationBoxplot(x)
# optimisationBoxplot(x, notch = TRUE)
