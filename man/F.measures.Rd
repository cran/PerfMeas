\name{F.measures}
\alias{F.measures}
\alias{F.measure.single}
\alias{F.measure.single.over.classes}
\alias{compute.mean.F.measure.single.over.classes}

\title{
F-measures
}
\description{
Set of functions to compute the F-measure, precision, recall, specificity and 0/1 loss accuracy.
}
\usage{
F.measure.single(pred, labels)
F.measure.single.over.classes(target, predicted, g, root = "00")
compute.mean.F.measure.single.over.classes(y)
}

\arguments{
   \item{pred}{
vector of the predicted labels. 0  stands for negative and 1 for positive
}
  \item{labels}{
vector of the true labels. 0  stands for negative and 1 for positive
}
  \item{target}{
matrix with the target multilabels. 0  stands for negative and 1 for positive.
Rows correspond to examples and columns to classes.
}
  \item{predicted}{
matrix with the predicted multilabels. 0  stands for negative and 1 for positive.
Rows correspond to examples and columns to classes.
}
  \item{g}{
graph of the classes (object of class graphNEL, package \pkg{graph}). 
If missing, no per level results are computed.
}
  \item{root}{
the name of the root node (def. "00") of the graph g.
}
  \item{y}{
a list of lists. The components of the outer list is a list returned from the function \code{F.measure.single.over.classes}
}
}
\details{
\code{F.measure.single} computes the F.score, precision, recall, specificity and accuracy for a single class.
  
  \code{F.measure.single.over.classes}  computes  precision, recall, specificity, accuracy and F-measure for a set of classes. In particualr it computes the correponding average values across classes, the average values across levels of the hierarchy of the classes (if any), and the values of the measures for each class. Note that if there is no hierarchy between classes (represented by the graph g), you can miss the g parameter and no per-level values are computed.
  
  \code{compute.mean.F.measure.single.over.classes} computes means across folds \cr 
  of F.measure.single.over.classes.   
  This function could be useful in cross-validated or multiple hold-out experimental settings.
}
\value{
  \code{F.measure.single} returns a named numeric vector with six elements:
     \item{- P }{precision}
	 \item{- R }{recall (sensitivity)}
	 \item{- S }{specificity}
	 \item{- F }{F measure}
	 \item{- A }{0/1 loss accuracy}
	 \item{- npos }{number of positive examples}
  
  \code{F.measure.single.over.classes}  returns  a list with three elements:
	\item{- average }{a named vector with the average precision, recall, specificity, F-measure, accuracy and average number of positive examples across classes.}
	\item{- per.level}{a named matrix with average  precision, recall, specificity, F-measure and accuracy for each level of the hierarchy. Named rows correspond to levels, named columns correspond respectively to precision, recall, specificity, F-measure, accuracy and number of positive examples.}
	\item{- per.class }{a named matrix with  precision, recall, specificity, F-measure, accuracy and number of positive examples for each class. Named rows correspond to classes, named columns correspond respectively to precision, recall, specificity, F-measure, accuracy and and number of positive examples.}
	
  \code{compute.mean.F.measure.single.over.classes}  returns a list obtained by averaging the results across folds of the input y. 
     The components are:
	\item{- average }{a named vector with the average precision, recall, specificity, F-measure and accuracy across classes across folds.}
	\item{- per.level}{a named matrix with average  precision, recall, specificity, F-measure and accuracy for each level of the hierarchy across folds. Named rows correspond to levels, named columns correspond respectively to precision, recall, specificity, F-measure and accuracy}
	\item{- per.class }{a named matrix with  precision, recall, specificity, F-measure and accuracy for each class across folds. Named rows correspond to classes, named columns correspond respectively to precision, recall, specificity, F-measure and accuracy.}
	
}

\seealso{
\code{\link{AUC.measures}}, \code{\link{PXR}}
}
\examples{
# preparing pseudo-random predictions and target-labels for examples: 100 examples
# and 10 classes
Scores <- matrix(runif(1000),nrow=100);
Targets <- Pred <- matrix(integer(1000),nrow=100);
Targets[Scores>0.5] <- 1;
# adding noise to scores
Scores <- Scores + matrix(rnorm(1000, sd=0.3),nrow=100);
Pred[Scores>0.5] <- 1;
colnames(Pred) <-colnames(Targets) <- LETTERS[1:10];
# getting predictions and labels of class "A"
pred <- Pred[,"A"];
labels <- Targets[,"A"];
# F.score and other metrics for a single class
F.measure.single(pred,labels);
# F.score and other metrics for the 10 classes
F.measure.single.over.classes(Targets, Pred);
}
