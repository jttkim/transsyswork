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


integratedPanelsEps(jdata.aggr, jGeneGroupOriginal, "jdata_aggr_%s_panel.eps");
integratedPanelsEps(jdata.aggr.mnorm, jGeneGroupOriginal, "jdata_aggr_mnorm_%s_panel.eps");
integratedPanelsEps(jdata.aggr.log, jGeneGroupOriginal, "jdata_aggr_log_%s_panel.eps");
integratedPanelsEps(jdata, jGeneGroupOriginal, "jdata_%s_panel.eps");
epsdevice("jdata_raw_groupdegree.eps");
par(cex.axis = 0.7);
groupDegreeBoxplot(jdata, jGeneGroupOriginal);
dev.off();
r <- unlist(getMostRepresentative(jdata, jGeneGroup));
cat(file = "mostrepresentative.txt", sprintf("%s: %s", names(r), r), sep = "\n");



# optimisationPlotEps("testjsmall%s_grad_squaresum_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_grad_squaresum_logratio_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_grad_correlation_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_grad_correlation_logratio_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_squaresum_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_squaresum_logratio_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_correlation_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_correlation_logratio_tfhypplus_r50");

optimisationPlotEps("testjfull%s_grad_squaresum_tfhypplus_r50");
optimisationPlotEps("testjfull%s_grad_squaresum_logratio_tfhypplus_r50");
optimisationPlotEps("testjfull%s_grad_correlation_tfhypplus_r50");
optimisationPlotEps("testjfull%s_grad_correlation_logratio_tfhypplus_r50");
optimisationPlotEps("testjfull%s_sastepadapt_squaresum_tfhypplus_r50");
optimisationPlotEps("testjfull%s_sastepadapt_squaresum_logratio_tfhypplus_r50");
optimisationPlotEps("testjfull%s_sastepadapt_correlation_tfhypplus_r50");
optimisationPlotEps("testjfull%s_sastepadapt_correlation_logratio_tfhypplus_r50");

profilePlotEps("testjfullprg_grad_correlation_logratio_tfhypplus_r50_profiles", unlist(jGeneGroup));
