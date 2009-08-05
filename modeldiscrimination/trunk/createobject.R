library(Biobase)


createobject <- function(expr, pheno)
{
  expressionSet <- readExpressionSet(exprsFile = expr, phenoDataFile = pheno);
  return(expressionSet);
}

