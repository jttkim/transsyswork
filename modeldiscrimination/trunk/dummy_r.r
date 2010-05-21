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
