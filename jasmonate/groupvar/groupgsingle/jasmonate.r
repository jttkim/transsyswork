### Read the coi1 / jasmonate data and organise it into a data frame.
### Column headers are expected in the format
###
###     <mutant>_<treatment>_t<time>_s<number>
### where
###
### * mutant is wt (wild type) or coi1
### * treatment is control (no treatment), jasmonate or wounding
### * time is the time from application of the treatment (always 0 for
###   controls, greater than 0 for treated plants
### * number is the number of the repetition of the experiment
###
### Attributes are provided to facilitate further processing and
### analysis. Each attribute has the one element to describe each
### row of the data frame, and the names of the attribute correspond
### to the colnames of the data frame.
###
### * mutant: "wt" or "coi1"
### * condition: "control", "wounding" or "jasmonate"
### * mtime: the measurement time in minutes


library(jtkstuff);


source("transopt.r");


hitReturn <- function(msg)
{
  invisible(readline(sprintf("%s -- hit return", msg)));
}


readJData <- function(filedescriptor)
{
  jdata <- read.table(filedescriptor, sep = ",", header = TRUE);
  jdata <- jdata[as.character(jdata[["affy_number"]]) != "", ];
  l <- strsplit(as.character(jdata[[2]]), " ");
  rownames(jdata) <- sapply(l, function(x) { return(x[1]); });
  jdata <- jdata[, 3:ncol(jdata)]
  mutant <- character(ncol(jdata));
  for (mname in c("wt", "coi1"))
  {
    r <- sprintf("%s_", mname);
    mutant[regexpr(r, colnames(jdata)) > -1] <- mname;
  }
  names(mutant) <- colnames(jdata);
  attr(jdata, "mutant") <- mutant;
  condition <- character(ncol(jdata));
  for (cname in c("control", "wounding", "jasmonate"))
  {
    r <- sprintf("_%s_", cname);
    condition[regexpr(r, colnames(jdata)) > -1] <- cname;
  }
  names(condition) <- colnames(jdata);
  attr(jdata, "condition") <- condition;
  # time of measurement
  mtime <- numeric(ncol(jdata));
  m <- regexpr("_t[0-9]+_", colnames(jdata));
  for (i in seq(along = m))
  {
    if (m[i] == -1)
    {
      mtime[i] <- NA;
    }
    else
    {
      j <- m[i] + 2;
      k <- m[i] + attr(m, "match.length")[i] - 2;
      mtime[i] <- as.numeric(substr(colnames(jdata)[i], j, k));
    }
  }
  names(mtime) <- colnames(jdata);
  attr(jdata, "mtime") <- mtime;
  return(jdata);
}


### Collapse columns that are identical in mutant, condition and time into
### one column containing the means.

aggregateData <- function(rawdata)
{
  # rawdata[rawdata < 0.0] <- NA;
  mutant <- attr(rawdata, "mutant");
  mtime <- attr(rawdata, "mtime");
  condition <- attr(rawdata, "condition");
  mutant.unique <- sort(unique(mutant));
  mtime.unique <- sort(unique(mtime));
  condition.unique <- sort(unique(condition));
  a <- data.frame(rep(NA, nrow(rawdata)), row.names = rownames(rawdata));
  i <- 1;
  a.mutant <- "???";
  a.mtime <- -1;
  a.condition <- "???"
  for (mut in mutant.unique)
  {
    for (cond in condition.unique)
    {
      for (mt in mtime.unique)
      {
        cols <- (mutant == mut) & (mtime == mt) & (condition == cond);
        if (sum(cols) > 0)
        {
          if (sum(cols) == 1)
          {
            a[, i] <- rawdata[, cols];
          }
          else
          {
            a[, i] <- rowMeans(rawdata[, cols], na.rm = TRUE);
          }
          colnames(a)[i] <- sprintf("%s_%s_t%d_mean", mut, cond, as.integer(mt));
          a.mutant[i] <- mut;
          a.condition[i] <- cond;
          a.mtime[i] <- mt;
          i <- i + 1;
        }
      }
    }
  }
  names(a.mutant) <- colnames(a);
  attr(a, "mutant") <- a.mutant;
  names(a.condition) <- colnames(a);
  attr(a, "condition") <- a.condition;
  names(a.mtime) <- colnames(a);
  attr(a, "mtime") <- a.mtime;
  return(a);
}


### Normalise by dividing each row by its mean.

normaliseMeans <- function(x)
{
  m <- rowMeans(jdata, na.rm = TRUE);
  for (i in 1:ncol(x))
  {
    x[[i]] <- x[[i]] / m;
  }
  return(x);
}


### Return a subset of columns, preserving the mutant,
### condition and mtime attributes.

arraySubset <- function(x, subrange)
{
  s <- x[, subrange, drop = FALSE];
  for (a in c("mutant", "condition", "mtime"))
  {
    attr(s, a) <- attr(x, a)[subrange];
  }
  return(s);
}


### Return a subset of genes (rows), preserving the
### mutant, condition and mtime attributes.

geneSubset <- function(x, geneList)
{
  s <- x[geneList, ];
  for (a in c("mutant", "condition", "mtime"))
  {
    attr(s, a) <- attr(x, a);
  }
  return(s);
}


### compute the means of gene groups.

groupMeans <- function(x, geneGroup)
{
  a <- data.frame(rep(NA, length(geneGroup)), row.names = names(geneGroup));
  for (i in 1:ncol(x))
  {
    a[[i]] <- rep(NA, length(geneGroup));
  }
  for (gg in names(geneGroup))
  {
    a[gg, ] <- colMeans(x[rownames(x) %in% geneGroup[[gg]], ]);
  }
  colnames(a) <- colnames(x);
  for (cn in c("mutant", "condition", "mtime"))
  {
    attr(a, cn) <- attr(x, cn);
  }
  return(a);
}


### Add a global offset to the (absolute) expression data.

addOffset <- function(data, o)
{
  odata <- data + o;
  for (a in c("mutant", "condition", "mtime"))
  {
    attr(odata, a) <- attr(data, a);
  }
  return(odata);
}


### Return the log ratios in a data frame.

logRatios <- function(data)
{
  mutant <- attr(data, "mutant");
  mutant.unique <- unique(mutant);
  condition <- attr(data, "condition");
  ctrl <- condition == "control";
  noncontrol <- condition != "control";
  lr <- data;
  for (m in mutant.unique)
  {
    mcols <- mutant == m;
    ccol <- data[, mcols & (condition == "control")];
    for (i in which(mcols))
    {
      lr[[i]] <- log2(lr[[i]] / ccol);
    }
  }
  attr(lr, "mutant") <- mutant;
  attr(lr, "condition") <- condition;
  attr(lr, "mtime") <- attr(data, "mtime");
  return(lr);
}


### Integrate the controls as measurements at t = 0 to the
### corresponding time series

integrateControls <- function(data)
{
  controlColumns <- attr(data, "condition") == "control";
  controlIndices <- which(controlColumns);
  otherColumns <- !controlColumns;
  otherIndices <- which(otherColumns);
  mutants.unique <- unique(attr(data, "mutant")[otherIndices]);
  l <- list();
  newMutant <- character();
  newCondition <- character();
  newMtime <- double();
  for (mut in mutants.unique)
  {
    mutControlColumns <- (attr(data, "mutant") == mut) & controlColumns;
    mutOtherColumns <- (attr(data, "mutant") == mut) & otherColumns;
    conditions.unique <- unique(attr(data, "condition")[mutOtherColumns]);
    for (cond in conditions.unique)
    {
      for (controlCol in colnames(data)[mutControlColumns])
      {
        newCol <- sub("control", cond, controlCol);
        l[[newCol]] <- data[[controlCol]];
        newMutant[newCol] <- mut;
        newCondition[newCol] <- cond;
        newMtime[newCol] <- 0;
      }
      condOtherColumns <- (attr(data, "condition") == cond) & mutOtherColumns;
      for (newCol in colnames(data)[condOtherColumns])
      {
        l[[newCol]] <- data[[newCol]];
        newMutant[newCol] <- mut;
        newCondition[newCol] <- cond;
        newMtime[newCol] <- attr(data, "mtime")[newCol];
      }
    }
  }
  d <- as.data.frame(l);
  rownames(d) <- rownames(data);
  attr(d, "mutant") <- newMutant;
  attr(d, "condition") <- newCondition;
  attr(d, "mtime") <- newMtime;
  return(d);
}


atopExpressionLabel <- function(rawLabel, splitString = "_")
{
  a <- expression();
  for (l in rawLabel)
  {
    s <- unlist(strsplit(l, splitString));
    if (length(s) > 4)
    {
      stop("too many components");
    }
    if (length(s) < 4)
    {
      s <- c(rep("", 4 - length(s)), s);
    }
    s1 <- s[1];
    s2 <- s[2];
    s3 <- s[3];
    s4 <- s[4];
    a <- c(a, as.expression(substitute(atop(atop(s1, s2), atop(s3, s4)))));
  }
  return(a);
}


stackedLabel <- function(rawLabel, splitString = "_")
{
  return(sub(splitString, "\n", rawLabel));
}


xaxisParameters <- function(data, xstyle)
{
  x.at <- NULL;
  x.labels <- TRUE;
  if (xstyle == "mtimes")
  {
    x <- attr(data, "mtime")[1:ncol(data)];
  }
  else
  {
    x <- 1:ncol(data);
    x.at <- x;
    if (xstyle == "mutant")
    {
      x.labels <- attr(data, "mutant");
    }
    else if (xstyle == "condition")
    {
      x.labels <- attr(data, "condition");
    }
    else if (xstyle == "stacked")
    {
      if (is.null(attr(data, "mutant")) || is.null(attr(data, "condition")) || is.null(attr(data, "mtime")))
      {
        x.labels <- stackedLabel(colnames(data));
      }
      else
      {
        conditionLabel <- attr(data, "condition");
        conditionLabel[conditionLabel == "control"] <- "ctrl";
        conditionLabel[conditionLabel == "jasmonate"] <- "ja";
        conditionLabel[conditionLabel == "wounding"] <- "wound";
        x.labels <- sprintf("%s\n%s\n%d", attr(data, "mutant"), conditionLabel, as.integer(attr(data, "mtime")));
      }      
    }
    else if (xstyle == "atop")
    {
      x.labels <- atopExpressionLabel(colnames(data));
    }
    else
    {
      if (is.null(attr(data, "mutant")) || is.null(attr(data, "condition")) || is.null(attr(data, "mtime")))
      {
        x.labels <- colnames(data);
      }
      else
      {
        conditionLabel <- attr(data, "condition");
        conditionLabel[conditionLabel == "control"] <- "c";
        conditionLabel[conditionLabel == "jasmonate"] <- "j";
        conditionLabel[conditionLabel == "wounding"] <- "w";
        x.labels <- sprintf("%s,%s,%d", attr(data, "mutant"), conditionLabel, as.integer(attr(data, "mtime")));
      }
    }
  }
  return(list(x = x, x.at = x.at, x.labels = x.labels));
}


### Display profiles by plotting one point for each
### measurement.

pointProfiles <- function(data, geneNames, x, curveLabels = FALSE, ...)
{
  n <- ncol(data);
  if (is.null(geneNames))
  {
    geneNames <- rownames(data);
  }
  # print(x);
  for (gene in geneNames)
  {
    y <- data[gene, 1:n];
    # print(sprintf("length(x) = %d, length(y) = %d", length(x), length(y)));
    points(x, y, ...);
    if (curveLabels)
    {
      text(x[n], y[n], gene, pos = 4, ...);
    }
  }
}


### Display profiles by plotting lines connecting measurements.

lineProfiles <- function(data, geneNames, x, curveLabels = FALSE, ...)
{
  n <- ncol(data);
  if (is.null(geneNames))
  {
    geneNames <- rownames(data);
  }
  for (gene in geneNames)
  {
    y <- data[gene, 1:n];
    lines(x, y, ...);
    if (curveLabels)
    {
      text(x[n], y[n], gene, pos = 4, offset = 0.0, ...);
    }
  }
}


plotProfiles <- function(data, geneNames = NULL, plotFunc = pointProfiles, xlim = NULL, ylim = NULL, xstyle = "mtimes", curveLabels = FALSE, labelSpace = curveLabels, main = NULL, ...)
{
  n <- ncol(data);
  xaxis.params <- xaxisParameters(data, xstyle);
  if (is.null(geneNames))
  {
    geneNames <- rownames(data);
  }
  if (is.null(xlim))
  {
    xlim <- c(min(xaxis.params$x), max(xaxis.params$x));
    if (xlim[2] == xlim[1])
    {
      xlim[2] <- xlim[2] + 1.0;
    }
    if (labelSpace)
    {
      ## FIXME: should use max. stringwidth of gene names here
      r <- xlim[2] - xlim[1];
      ## 1.3 is a "magic" value, intended to reflect the ratio of total width to plotarea width
      w <- max(strwidth(rownames(data), "figure")) * 1.3;
      lw <- w * r / (1.0 - w);
      xlim[2] <- xlim[2] + lw;
    }
  }
  if (is.null(ylim))
  {
    ylim <- range(data[geneNames, 1:n]);
  }
  ## print(xlim);
  ## print(ylim);
  ## print(geneNames);
  plot.new();
  plot.window(xlim, ylim, ...);
  plotFunc(data, geneNames, xaxis.params$x, curveLabels, ...);
  axis(1, at = xaxis.params$x.at, labels = xaxis.params$x.labels, ...);
  axis(2, ...);
  box(...);
  title(main = main, ...);
  return(xaxis.params);
}


plotMutants <- function(data, geneNames = NULL, plotFunc = pointProfiles, xstyle = "mtimes", ...)
{
  ## palette(c("red", "green", "blue", "cyan", "magenta", "black", "yellow"));
  if (is.null(geneNames))
  {
    geneNames <- rownames(data);
  }
  mutant <- unique(attr(data, "mutant"));
  palette(rainbow(length(mutant)));
  xaxis.params <- plotProfiles(data, geneNames, plotFunc = plotFunc, xstyle = xstyle, ...);
  for (i in seq(along = mutant))
  {
    subdata <- arraySubset(data, attr(data, "mutant") == mutant[i]);
    xp <- xaxisParameters(subdata, xstyle);
    plotFunc(subdata, geneNames, xp$x, col = i, ...);
  }
}


plotAllGroups <- function(data, geneGroupList, xstyle = "mtimes", waitFunc = hitReturn)
{
  for (g in names(geneGroupList))
  {
    plotMutants(data, geneGroupList[[g]], xstyle = xstyle);
    waitFunc(g);
  }
}


plotAllProfiles <- function(data, geneGroupList, plotFunc = pointProfiles, xstyle = "x", ...)
{
  ## palette(c("red", "green", "blue", "cyan", "magenta", "black", "yellow"));
  palette(rainbow(length(geneGroupList)));
  xaxis.params <- plotProfiles(data, rownames(data), plotFunc, xstyle = xstyle, ...);
  for (i in seq(along = geneGroupList))
  {
    plotFunc(data, geneGroupList[[i]], xaxis.params$x, col = i, ...);
  }
}


plotHighlightProfiles <- function(data, geneGroupList, plotFunc = pointProfiles, xstyle = "mtimes", highlightColour = "red", waitFunc = hitReturn, ...)
{
  for (g in names(geneGroupList))
  {
    xaxis.params <- plotProfiles(data, NULL, plotFunc, xstyle = xstyle, ...);
    plotFunc(data, geneGroupList[[g]], xaxis.params$x, col = highlightColour, ...);
    waitFunc(g);
  }
}


plotHighlightPanels <- function(data, geneGroupList, plotFunc = pointProfiles, xstyle = "atop", highlightColour = "red", waitFunc = hitReturn, ...)
{
  mutant <- attr(data, "mutant");
  mutant.unique <- unique(mutant);
  condition = attr(data, "condition");
  condition.unique <- unique(condition);
  opar <- par(no.readonly = TRUE);
  par(mfcol = c(length(condition.unique), length(mutant.unique)));
  for (g in names(geneGroupList))
  {
    for (mut in mutant.unique)
    {
      subdata.mut <- arraySubset(data, mutant == mut);
      for (cond in condition.unique)
      {
        subdata.mut.cond <- arraySubset(subdata.mut, attr(subdata.mut, "condition") == cond)
        xaxis.params <- plotProfiles(subdata.mut.cond, NULL, plotFunc = plotFunc, xstyle = xstyle, main = sprintf("group %s, %s, %s", g, mut, cond), curveLabels = FALSE, labelSpace = TRUE, ...);
        plotFunc(subdata.mut.cond, geneGroupList[[g]], xaxis.params$x, col = highlightColour, curveLabels = TRUE, ...);
      }
    }
    waitFunc(g);
  }
  par(opar);
}


getMeans <- function(x)
{
  mutants <- unique(attr(x, "mutant"));
  conditions <- unique(attr(x, "condition"));
  l <- list();
  l[["all"]] <- rowMeans(x);
  for (mut in mutants)
  {
    l[[mut]] <- rowMeans(x[, attr(x, "mutant") == mut]);
  }
  for (cond in conditions)
  {
    l[[cond]] <- rowMeans(x[, attr(x, "condition") == cond]);
  }
  for (mut in mutants)
  {
    for (cond in conditions)
    {
      l[[sprintf("%s.%s")]] <- rowMeans(x[, attr(x, "mutant") == mut & attr(x, "condition") == cond])
    }
  }
  return(as.data.frame(l));
}


getDifferential <- function(x, cols)
{
  xplus <- x[, cols];
  xminus <- x[, !cols];
  mplus <- rowMeans(xplus);
  sdplus <- apply(xplus, 1, sd);
  mminus <- rowMeans(xminus);
  sdminus <- apply(xminus, 1, sd);
  return((mplus - mminus) / (sdplus + sdminus));
}


writeJasmonateObjectiveFile <- function(filename, data)
{
  file <- file(filename);
  open(file, open = "w");
  cat(sprintf("mutant: %s\n", paste(attr(data, "mutant"), sep = "", collapse = " ")), file = file);
  cat(sprintf("condition: %s\n", paste(attr(data, "condition"), sep = "", collapse = " ")), file = file);
  v <- sprintf("%7.2f", attr(data, "mtime"));
  cat(sprintf("mtime: %s\n", paste(v, sep = "", collapse = " ")), file = file);
  for (geneName in rownames(data))
  {
    s <- sprintf("%s:", geneName);
    v <- sprintf("%1.17e", data[geneName, ]);
    l <- paste(c(s, v), sep = "", collapse = " ");
    cat(sprintf("%s\n", l), file = file);
  }
  close(file);
}


groupDegrees <- function(x)
{
  l <- list();
  l[["a"]] <- getDifferential(x, attr(x, "mutant") == "wt" & attr(x, "condition") != "control");
  l[["b"]] <- getDifferential(x, attr(x, "condition") != "control");
  l[["bc"]] <- -getDifferential(x, attr(x, "mutant") == "wt" & attr(x, "condition") != "control");
  l[["c"]] <- -getDifferential(x, attr(x, "mutant") == "wt" & attr(x, "condition") != "control");
  l[["d"]] <- getDifferential(x, attr(x, "condition") == "wounding");
  l[["e"]] <- -getDifferential(x, attr(x, "condition") != "control");
  l[["f"]] <- getDifferential(x, attr(x, "mutant") == "wt" & attr(x, "condition") == "jasmonate");
  l[["g"]] <- -getDifferential(x, attr(x, "mutant") == "wt" & attr(x, "condition") == "wounding");
  return(as.data.frame(l));
}


getProfileList <- function(profilePairs, geneName, profileSource)
{
  sourcecol <- which(colnames(profilePairs) == "source");
  if (length(sourcecol) != 1)
  {
    stop("source column not found");
  }
  profileColnames <- colnames(profilePairs)[(sourcecol + 1):ncol(profilePairs)];
  l <- list();
  for (i in which((profilePairs[["gene"]] == geneName) & (profilePairs[["source"]] == profileSource)))
  {
    restartIndex <- profilePairs[i, "restart_index"];
    p <- as.numeric(profilePairs[i, (sourcecol + 1):ncol(profilePairs)]);
    names(p) <- profileColnames;
    l[[as.character(restartIndex)]] <- p;
  }
  return(l);
}


groupDegreeBoxplot <- function(x, geneGroup, ...)
{
  d <- groupDegrees(x);
  l <- list();
  groupScore <- 0.0;
  for (groupName in names(geneGroup))
  {
    g <- rownames(d) %in% geneGroup[[groupName]];
    ## print(g);
    gplus <- d[g, groupName];
    gminus <- d[!g, groupName];
    l[[sprintf("%s+", groupName)]] <- gplus;
    l[[sprintf("%s-", groupName)]] <- gminus;
    if ((length(gplus) > 1) && (length(gminus) > 1))
    {
      groupScore <- groupScore + (mean(gplus) - mean(gminus)) / (sd(gplus) + sd(gminus));
    }
  }
  ## print(l);
  boxplot(l, ...);
  return(groupScore);
}


getMostRepresentative <- function(x, geneGroup)
{
  d <- groupDegrees(x);
  l <- list();
  for (groupName in names(geneGroup))
  {
    g <- rownames(d) %in% geneGroup[[groupName]];
    m <- max(d[g, groupName]);
    i <- which(d[[groupName]] == m);
    l[[groupName]] <- rownames(d)[i];
    message(sprintf("%s: %s (%f)", groupName, rownames(d)[i], m));
  }
  return(l);
}


getGenesWithProfileType <- function(profileData, profileType)
{
  return(as.character(unique(profileData[profileData[["source"]] == profileType, "gene"])));
}


getComparableGenes <- function(profileData)
{
  return(intersect(getGenesWithProfileType(profileData, "empirical"), getGenesWithProfileType(profileData, "simulated")));
}


# plot sets of profiles, empirical in red and simulated in green

plotProfileSets <- function(profilePairs, geneNames = NULL, main = NULL, waitFunc = hitReturn, ...)
{
  if (is.null(geneNames))
  {
    geneNames <- unique(as.character(profilePairs[["gene"]]));
  }
  for (geneName in geneNames)
  {
    elist <- getProfileList(profilePairs, geneName, "empirical");
    slist <- getProfileList(profilePairs, geneName, "simulated");
    y <- c(unlist(elist), unlist(slist));
    ylim <- c(min(y), max(y));
    if (length(elist) > 0)
    {
      y <- elist[[1]];
    }
    else if (length(slist) > 0)
    {
      y <- slist[[1]];
    }
    else
    {
      y <- NULL;
    }
    if (is.null(y))
    {
      message(sprintf("no profiles found for gene %s", geneName));
    }
    else
    {
      x <- 0:(length(y) - 1);
      if (is.null(main))
      {
        m = geneName;
      }
      else
      {
        m = main;
      }
      plot(x, y, ylim = ylim, axes = FALSE, type = "l", main = m);
      axis(1, at = x, labels = atopExpressionLabel(names(y)), ...);
      axis(2, ...);
      box(...);
      for (s in slist)
      {
        lines(x, s, col = "green");
      }
      for (e in elist)
      {
        lines(x, e, col = "red");
      }
    }
    waitFunc(geneName);
  }
}


plotProfileSetsCorr <- function(profilePairs, geneNames = NULL, waitFunc = hitReturn)
{
  if (is.null(geneNames))
  {
    geneNames <- unique(as.character(profilePairs[["gene"]]));
  }
  for (geneName in geneNames)
  {
    elist <- getProfileList(profilePairs, geneName, "empirical");
    slist <- getProfileList(profilePairs, geneName, "simulated");
    if (!all(names(elist) == names(slist)))
    {
      message(sprintf("cannot associate empirical and simulated profiles for gene %s", geneName));
    }
    else
    {
      x <- unlist(elist);
      y <- unlist(slist);
      plot(x, y, type = "p", xlim = c(min(x), max(x)), ylim = c(min(y), max(y)), main = geneName, xlab = "empirical", ylab = "simulated");
    }
    waitFunc(geneName);
  }
}


readOptimisationLogs <- function(fnameFormat, fileOpenFunc = NULL)
{
  return(readFrames(c("prg", "blank", sprintf("rnd%03d", 1:20)), fnameFormat, fileOpenFunc));
}


plotOptimisationLogs <- function(optLog, colNames = c("obj", "stepsize"), waitfunc = hitReturn)
{
  for (dfname in names(optLog))
  {
    for (colName in colNames)
    {
      main = sprintf("%s, %s", dfname, colName);
      plotListRainbow(tapply(optLog[[dfname]][[colName]], optLog[[dfname]][["restart_index"]], function(x) {x}), main = main);
      waitFunc(main);
    }
  }
}


optimisationBoxplot <- function(optLog, ...)
{
  finalObjList <- lapply(optLog, function(x) {getFinalStats(x)[["obj"]];});
  for (n in names(finalObjList))
  {
    print(sprintf("%s (%d samples): %f +/- %f", n, length(finalObjList[[n]]), mean(finalObjList[[n]]), sd(finalObjList[[n]])));
  }
  boxplot(finalObjList, ...);
}


generateJdata <- function()
{
  jdata <- readJData("coi1_affy_mod_repaired.csv");
  jGeneGroupOriginal <- list(
                             a = c("At1g18710", "At1g23080", "At1g74020", "At2g34810", "At4g11290", "At4g22880", "At4g23600", "At5g13930", "At5g19890", "At5g57090"),
                             b = c("At1g17740", "At2g21130", "At2g28900", "At2g42530", "At4g24360", "At5g63790", "At4g35770"),
                             c = c("At3g23050", "At4g35770"),
                             d = c("At2g23320", "At2g29500", "At2g38470", "At4g11280"),
                             e = c("At1g20440", "At2g26690", "At2g40900", "At3g61890", "At4g34000", "At4g37760", "At5g06760"),
                             f = c("At1g49570", "At4g11320"),
                             g = c("At4g17880", "At5g62470")
                             );
  jGeneGroupBc <- list(
                       a = c("At1g18710", "At1g23080", "At1g74020", "At2g34810", "At4g11290", "At4g22880", "At4g23600", "At5g13930", "At5g19890", "At5g57090"),
                       b = c("At1g17740", "At2g21130", "At2g28900", "At2g42530", "At4g24360", "At5g63790"),
                       bc = "At4g35770",
                       c = "At3g23050",
                       d = c("At2g23320", "At2g29500", "At2g38470", "At4g11280"),
                       e = c("At1g20440", "At2g26690", "At2g40900", "At3g61890", "At4g34000", "At4g37760", "At5g06760"),
                       f = c("At1g49570", "At4g11320"),
                       g = c("At4g17880", "At5g62470")
                       );
  jGeneGroup <- list(
                     a = c("At1g18710", "At1g23080", "At1g74020", "At2g34810", "At4g11290", "At4g22880", "At4g23600", "At5g13930", "At5g19890", "At5g57090"),
                     f = c("At1g49570", "At4g11320"),
                     c = c("At3g23050", "At4g35770"),
                     g = c("At4g17880", "At5g62470"),
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

                                        # writeJasmonateObjectiveFile("jdata_aggr.txt", jdata.aggr);
  writeJasmonateObjectiveFile("jdata_aggr_repaired.txt", jdata.aggr);
}

jGeneGroup <- list(
                   a = c("At1g18710", "At1g23080", "At1g74020", "At2g34810", "At4g11290", "At4g22880", "At4g23600", "At5g13930", "At5g19890", "At5g57090"),
                   f = c("At1g49570", "At4g11320"),
                   c = c("At3g23050", "At4g35770"),
                   g = c("At4g17880", "At5g62470"),
                   b = c("At1g17740", "At2g21130", "At2g28900", "At2g42530", "At4g24360", "At5g63790"),
                   e = c("At1g20440", "At2g26690", "At2g40900", "At3g61890", "At4g34000", "At4g37760", "At5g06760"),
                   d = c("At2g23320", "At2g29500", "At2g38470", "At4g11280")
                   );
