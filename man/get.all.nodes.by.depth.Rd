\name{get.all.nodes.by.depth}
\alias{get.all.nodes.by.depth}

\title{
Getting nodes by their depth
}
\description{
Grouping classes by level in a given hierarchy. 
}
\usage{
get.all.nodes.by.depth(g, root = "00")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
graph of the classes (object of class graphNEL, package \pkg{graph}). 
}
  \item{root}{
name of the root node (def. 00)
}
}

\details{The minimum paths between the \dQuote{root} and all the other classes/nodes are computed.  Levels are numbered from 1 in increasing order by their distance from the \dQuote{root} class.}

\value{
a  list of the nodes, grouped w.r.t. the distance from the root.
The first element of the list corresponds to the nodes at distance 1, 
the second to nodes at distance 2 and so on.
}

