\name{graphics}
\alias{precision.recall.curves.plot}
\alias{performance.curves.plot}

\title{
Graphics function to plot precision/recall or f.score/recall curves
}
\description{
Function to plot multiple precision/recall or f.score/recall curves
}
\usage{
precision.recall.curves.plot(y, range=seq(from=0, to=1, by=0.1), 
curve.names=1:length(y), cex.val=0.6, f="", height=9, width=11,
col=c("black","red1","blue1","green1","darkgrey","brown1","yellow1","orange1",
"red4","blue4","green4","lightgrey","brown4","yellow4","orange4"),
line.type=1, leg=TRUE, pos=c(range[length(range)-2], range[length(range)]), do.grid=TRUE,
plot.precision=TRUE, trap.rule=TRUE)

performance.curves.plot(m, x.range=seq(from=0.1, to=1, by=0.1), 
y.range=c(0,1), curve.names=1:nrow(m), cex.val=0.6, f="", height=9, width=11,
col=c("black","red1","blue1","green1","darkgrey","brown1","yellow1","orange1",
"red4","blue4","green4","lightgrey","brown4","yellow4","orange4"), line.type=1, 
patch.type=1:16, leg=TRUE, pos=c(x.range[length(x.range)-2], y.range[2]), do.grid=TRUE,
x.label="Recall", y.label="Precision")
}

\arguments{
  \item{y}{
a list of lists. 
Each component list is a list returned from \code{precision.at.all.recall.levels}
that reports precision and recall results at different levels for different methods or tasks
}
  \item{range}{
numeric vector of the precision/recall values to be represented (def: values between 0 and 1 step 0.1)
}
  \item{curve.names}{
names of the compared methods to be reported in the legenda (def: numbers)
}
  \item{cex.val}{
magnification value for characters (def. 0.6)
}
  \item{f}{
file name. If is given, an encapsulated postscript file is created, otherwise the output is rendered on
       a window.
}
  \item{height}{
relative heigth of the graph (def. 9)
}
  \item{width}{
relative width of the graph (def. 11)
}
  \item{col}{
colors of the lines. 14 different colors are given as default, but any vector of color from colors() (package \pkg{graphics}) can be used. Colors are recycled if length(col) < length(y).
}
  \item{line.type}{
type of the line. Any valid vector of integer can be assigned (values between 1  and 6, see lty in par(), package \pkg{graphics} for details). Values are recycled if length(line.type) < length(y). Def.: 1 (solid lines).
}
\item{leg}{ boolean: if TRUE (def.) a legend is depicted.}
\item{pos}{coordinates of the position of the legend.}
\item{plot.precision}{ boolean: if TRUE (def.) precision/recall curves are plotted, othewise f-score/recall curves.}
\item{trap.rule}{ boolean: if TRUE (def.) the integral of the curves are computed.}
\item{m}{ a numeric matrix. Rows correspond to different methods and columns to precision or f-score at given recall values }
\item{x.range}{vector of the recall values to be represented}
\item{y.range}{vector with 2 elements: range of the precision/f-score to be represented }
\item{patch.type}{ numeric vector corresponding to the symbols to be plotted for different recall values (def. 1:16)}
\item{do.grid}{a bollen value indicating whether a grid will be plotted or not (def: TRUE)}
\item{x.label}{label of the abscissa (def: Recall)}
\item{y.label}{label of the ordinate (def: Precision)}
}

\details{
The function \code{precision.recall.curves.plot} plots multiple precision/recall curves (or f-score/recall curves, if the argument \code{plot.precision=FALSE}) by using the output of the function \code{precision.at.all.recall.levels}, that compute several precison/recall pairs by moving the threshold from the lowest to the highest score achieved by each example.

The function \code{performance.curves.plot} plots precision of F-score/recall curves ofr a predefined set of recall levels.
Thsi function can be useful to plot and to compare the average results between methods across multiple classes.

The curves can differ by color, type of line and for \code{performance.curves.plot} for each recall value a symbol can be also plotted.  A legend can be automatically constructed.
}
\value{
The functions output  a graphic file either on a window or on an encapsulated postscript file.
The function \code{precision.recall.curves.plot}, if the parameter \code{trap.rule = TRUE} (default), outputs a vector with the AUPRC (or the Area Under the F-score Curve if the parameter \code{plot.precision=FALSE}) for each curve.
}

\examples{
# loading  an example matrix of scores an the correponding table of classes
data(T);
data(Scores);
res=list();
classes=c(1,2,7,8)
# computing precison recall values
for (j in classes) res=c(res, list(precision.at.all.recall.levels(Scores[,j], T[,j])));
# plotting precision/recall curves
precision.recall.curves.plot(res,curve.names=colnames(T)[classes], 
pos=c(0.7,1), plot.precision=TRUE);
# black and white version
precision.recall.curves.plot(res,curve.names=colnames(T)[classes], pos=c(0.7,1), 
plot.precision=TRUE, line.type=1:4, col=1);
# plotting f-score/recall curves
precision.recall.curves.plot(res,curve.names=colnames(T)[classes], pos=c(0.7,1), 
plot.precision=FALSE);
# black and white version
precision.recall.curves.plot(res,curve.names=colnames(T)[classes], pos=c(0.7,1), 
plot.precision=TRUE, line.type=1:4, col=1);
}
