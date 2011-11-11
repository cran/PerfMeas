\name{AUPRC}
\alias{AUPRC}
\alias{trap.rule.integral}

\title{
Area Under the Precision Recall Curve
}
\description{
Functions to compute the Area Under the Precision Recall Curve (AUPRC) and the Area Under the F-score Recall Curve (AUFRC)
}
\usage{
AUPRC(z, comp.precision=TRUE)
trap.rule.integral(x,y)
}

\arguments{
 \item{z}{
a list of lists. The components of the outer list is a list returned from the function \code{precision.at.all.recall.levels} that reports precision, recall and f-score results at different levels for different methods or tasks.
}
  \item{comp.precision}{
boolean. It TRUE (default) the AUPRC is computed otherwise the area under the F-score curve is computed
}
  \item{x}{
vector of the x values in increasing order
}
  \item{y}{
vector of the corresponding y=f(x) values
}


}
\details{
\code{AUPRC} computes the Area Under the Precision Recall Curve or the Area Under the F-score Recall Curve (AUFRC) 
for multiple curves by using the output of the function \code{precision.at.all.recall.levels}.
 
The function  \code{trap.rule.integral} implements the trapezoidal rule of integration and can be used to compute the integral of any empirical function expressed as a set of pair values (a vector of \code{x} values and a vector of \code{y = f(x)} values). In particular if \code{x} is the recall (with values in ascending order) and \code{y} the corresponding precision, \code{trap.rule.integral} copmutes the AUPRC.
}
\value{
  \code{AUPRC} returns the value of the AUPRC (if the argument \code{comp.precision = TRUE}), otherwise the value of the AUFRC.
  
\code{trap.rule.integral}  returns the value of the integral.
	
}

\seealso{
\code{\link{AUC.measures}}, \code{\link{PXR}}
}
\examples{
# loading matrices of scores an correponding table of classes
data(T);
data(Scores);
res=list();
classes=1:10
# computing precision recall values
for (j in classes) res=c(res, list(precision.at.all.recall.levels(Scores[,j], T[,j])));
names(res)<-seq(0.1, 1, by=0.1);
# computing AUPRC
AUPRC (res, comp.precision=TRUE);
# computing AU F-score recall curve
AUPRC (res, comp.precision=TRUE);

# Loading precision at given recall levels for different methods
data(PrecRec);
# computing AUPRC for different methods
x <- seq(0.1, 1, by=0.1);
res <- numeric(nrow(PrecRec));
names(res) <- rownames(PrecRec);
for (i in 1:nrow(PrecRec))  
  res[i] <- trap.rule.integral(x, PrecRec[i,]);
print(res);

}

