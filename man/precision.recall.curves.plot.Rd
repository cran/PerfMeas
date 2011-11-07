\name{precision.recall.curves.plot}
\alias{precision.recall.curves.plot}

\title{
Plotting precisions at given recall levels
}
\description{
Function to plot multiple precision recall curves
}
\usage{
precision.recall.curves.plot(y, curve.names = 1:length(y), f = "", 
range = seq(from = 0, to = 1, by = 0.1), cex.val = 0.6, height = 9, 
width = 11, col = c("black", "red1", "blue1", "green1", "darkgrey", 
"brown1", "yellow1", "orange1", "red4", "blue4", "green4", 
"lightgrey", "brown4", "yellow4", "orange4"), line.type = 1, 
leg=TRUE, pos=c(range[length(range)-2], y=range[length(range)]))
}

\arguments{
  \item{y}{
a list of lists. 
Each component list is a list returned from \code{precision.at.all.recall.levels}
that reports precision and recall results at different levels for different methods or tasks
}
  \item{curve.names}{
names of the compared methods to be reported in the legenda (def: numbers)
}
  \item{f}{
file name. If is given, an encapsulated postscript file is created, otherwise the output is rendered on
       a window.
}
  \item{range}{
numeric vector of the precision/recall values to be represented (def: values between 0 and 1 step 0.1)
}
  \item{cex.val}{
magnification value for characters (def. 0.6)
}
  \item{height}{
realstive heigth of the graph (def. 9)
}
  \item{width}{
relstive width of the graph (def. 11)
}
  \item{col}{
colors of the lines. 14 different colors are given as default, but any vector of color from colors() (package \pkg{graphics})
 can be used. Colors are recycled if length(col) < length(y).
}
  \item{line.type}{
type of the line. Any valid vector of integer can be assigned (values between 1  and 6, see lty in par(), package \pkg{graphics} for details). Values are recycled if length(line.type) < length(y). Def.: 1 (solid lines).
}
\item{leg}{ boolean: if TRUE (def.) a legend is depicted.}
\item{pos}{coordinates of the position of the legend.}
}

\details{
It plot multiple curves with a given recall in abscissa and the corresponding precision in ordinate. The curves can differ by color and tpe of line. A legend is automatically constructed.
}
\value{
It outputs a graphics file either on a window or on an encapusleted postscript file.
}

\examples{
# preparing pseudo-random predictions and target-labels to plot results of 5 methods
n=500;
Scores <- matrix(runif(n),nrow=100);
Targets <- matrix(integer(n),nrow=100);
Targets[Scores>0.5] <- 1;
# adding noise to scores
Scores <- Scores + matrix(rnorm(n, sd=0.7),nrow=100);
colnames(Scores) <-colnames(Targets) <- paste("Method", 1:5);
res=list();
for (j in 1:(n/100)) 
   res=c(res, list(precision.at.all.recall.levels(Scores[,j], Targets[,j])));
precision.recall.curves.plot(res,curve.names=paste("Method", 1:5), pos=c(0.7,0.3));
}
