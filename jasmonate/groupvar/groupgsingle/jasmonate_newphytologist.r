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


profilePlotEps <- function(baseName, geneGroupList, ...)
{
  profileFile <- sprintf("%s.dat", baseName);
  pp <- read.table(profileFile, header = TRUE);
  for (geneGroupName in names(geneGroupList))
  {
    for (geneName in geneGroupList[[geneGroupName]])
    {
      epsFileName <- sprintf("%s_%s_%s.eps", baseName, geneGroupName, geneName);
      epsdevice(epsFileName);
      plotProfileSets(pp, geneName, main = sprintf("%s (group %s)", geneName, geneGroupName), waitFunc = function(x) { print(x); });
      dev.off();
    }
  }
}


optimisationPlotEps <- function(baseFormat)
{
  fnameFormat <- sprintf("%s_log.dat.gz", baseFormat);
  print(fnameFormat);
  l <- readOptimisationLogs(fnameFormat, fileOpenFunc = gzfile);
  epsFileName <- sprintf("newphytologist_%s_optbox.eps", sprintf(baseFormat, ""));
  epsdevice(epsFileName);
  par(cex.axis = 1.0);
  optimisationBoxplot(l, notch = TRUE);
  dev.off();
}


optimisationPlotEpsPolishedRandomRewired <- function(baseFormat, ylab)
{
  fnameFormat <- sprintf("%s_log.dat.gz", baseFormat);
  print(fnameFormat);
  l <- readOptimisationLogs(fnameFormat, fileOpenFunc = gzfile);
  epsFileName <- sprintf("newphytologist_%s_optbox.eps", sprintf(baseFormat, ""));
  epsdevice(epsFileName);
  par(cex.axis = 0.7, cex.lab = 1.4, mar = c(5, 5, 1, 1) + 0.1);
  optimisationBoxplot(l, notch = TRUE, names = c("p", "n", rep("r", 20)), xlab = "model", ylab = ylab);
  dev.off();
}


jdata <- readJData("coi1_affy_mod_final.csv");
## jGeneGroupOriginal <- list(
##                            a = c("At1g18710", "At1g23080", "At1g74020", "At2g34810", "At4g11290", "At4g22880", "At4g23600", "At5g13930", "At5g19890", "At5g57090"),
##                            b = c("At1g17740", "At2g21130", "At2g28900", "At2g42530", "At4g24360", "At5g63790", "At4g35770"),
##                            c = c("At3g23050", "At4g35770"),
##                            d = c("At2g23320", "At2g29500", "At2g38470", "At4g11280"),
##                            e = c("At1g20440", "At2g26690", "At2g40900", "At3g61890", "At4g34000", "At4g37760", "At5g06760"),
##                            f = c("At1g49570", "At4g11320"),
##                            g = c("At4g17880", "At5g62470")
##                            );
## jGeneGroupBc <- list(
##                      a = c("At1g18710", "At1g23080", "At1g74020", "At2g34810", "At4g11290", "At4g22880", "At4g23600", "At5g13930", "At5g19890", "At5g57090"),
##                      b = c("At1g17740", "At2g21130", "At2g28900", "At2g42530", "At4g24360", "At5g63790"),
##                      bc = "At4g35770",
##                      c = "At3g23050",
##                      d = c("At2g23320", "At2g29500", "At2g38470", "At4g11280"),
##                      e = c("At1g20440", "At2g26690", "At2g40900", "At3g61890", "At4g34000", "At4g37760", "At5g06760"),
##                      f = c("At1g49570", "At4g11320"),
##                      g = c("At4g17880", "At5g62470")
##                      );
jGeneGroup <- list(
                   a = c("At1g18710", "At1g23080", "At1g74020", "At2g34810", "At4g11290", "At4g22880", "At4g23600", "At5g13930", "At5g19890", "At5g57090"),
                   f = c("At1g49570", "At4g11320"),
                   c = c("At3g23050", "At4g35770", "At5g62470"),
                   g = c("At4g17880"),
                   b = c("At1g17740", "At2g21130", "At2g28900", "At2g42530", "At4g24360", "At5g63790"),
                   e = c("At1g20440", "At2g26690", "At2g40900", "At3g61890", "At4g34000", "At4g37760", "At5g06760"),
                   d = c("At2g23320", "At2g29500", "At2g38470", "At4g11280")
                   );
## a: genes that are induced by MeJA via COI1, and also by wounding, mediated by JA
## b: genes induced by MeJA and wounding, independently of COI1
## bc: "intriguingly, At4g35770 is also repressed by wounding and MeJa, in a COI1-dependent manner"
## c: repressed by wounding and MeJA in a COI1-dependent manner
## d: genes induced by wounding but not MeJA in a COI-independent manner
## e: MeJA and wound-mediated COI1-independent repression
## f: COI1-mediated, MeJA but not wound- induction
## g: repressed by wounding only in a COI-dependent manner

jdata.wt <- arraySubset(jdata, attr(jdata, "mutant") == "wt");
jdata.coi1 <- arraySubset(jdata, attr(jdata, "mutant") == "coi1");
jdata.mnorm <- normaliseMeans(jdata);
jdata.aggr <- aggregateData(jdata);
jdata.aggr.mnorm <- normaliseMeans(jdata.aggr);
jdata.aggr.log <- logRatios(addOffset(jdata.aggr, -2.0 * min(jdata.aggr)));
jdata.aggr.log.gmeans <- groupMeans(jdata.aggr.log, jGeneGroup);
jdata.gmeans <- groupMeans(jdata, jGeneGroup);
jdata.mnorm.gmeans <- groupMeans(jdata.mnorm, jGeneGroup);

integratedPanelsEps(jdata.aggr, jGeneGroupOriginal, "jdata_aggr_%s_panel.eps");
integratedPanelsEps(jdata.aggr.mnorm, jGeneGroupOriginal, "jdata_aggr_mnorm_%s_panel.eps");
integratedPanelsEps(jdata.aggr.log, jGeneGroupOriginal, "jdata_aggr_log_%s_panel.eps");
integratedPanelsEps(jdata, jGeneGroupOriginal, "jdata_%s_panel.eps");
## epsdevice("jdata_raw_groupdegree.eps");
## par(cex.axis = 0.7);
## groupDegreeBoxplot(jdata, jGeneGroupOriginal);
## dev.off();
## r <- unlist(getMostRepresentative(jdata, jGeneGroup));
## cat(file = "mostrepresentative.txt", sprintf("%s: %s", names(r), r), sep = "\n");

# testjfullrnd019_grad_tfhypplus_correlation_logratio_log.dat.gz

optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_tfhypplus_squaresum", expression(D[ssq]));
optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_tfhypplus_squaresum_logratio", expression(D[ssq]^"(l)"));
optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_tfhypplus_correlation", expression(D[corr]^"(l)"));
optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_tfhypplus_correlation_logratio", expression(D[corr]^"(l)"));

## optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_squaresum_tfhypplus_r50", expression(D[ssq]));
## optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_squaresum_logratio_tfhypplus_r50", expression(D[ssq]^"(l)"));
## optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_correlation_tfhypplus_r50", expression(D[corr]^"(l)"));
## optimisationPlotEpsPolishedRandomRewired("testjfull%s_grad_correlation_logratio_tfhypplus_r50", expression(D[corr]^"(l)"));


# optimisationPlotEps("testjsmall%s_grad_squaresum_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_grad_squaresum_logratio_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_grad_correlation_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_grad_correlation_logratio_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_squaresum_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_squaresum_logratio_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_correlation_tfhypplus_r50");
# optimisationPlotEps("testjsmall%s_sastepadapt_correlation_logratio_tfhypplus_r50");

## optimisationPlotEps("testjfull%s_sastepadapt_squaresum_tfhypplus_r50");
## optimisationPlotEps("testjfull%s_sastepadapt_squaresum_logratio_tfhypplus_r50");
## optimisationPlotEps("testjfull%s_sastepadapt_correlation_tfhypplus_r50");
## optimisationPlotEps("testjfull%s_sastepadapt_correlation_logratio_tfhypplus_r50");

profilePlotEps("testjfullprg_grad_tfhypplus_correlation_logratio_profiles", jGeneGroup);
