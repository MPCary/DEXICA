#' Data: ANN model for use in component partitioning
#'
#' An artificial neural network (ANN) trained using simulated data to find the
#' optimal thresholds for partitioning components into discrete sets.
#'
#' @format An object of class "nn" built using the neuralnet package
#'
#' @source See PMID:(TBD) for details of ann construction.
"ann.partition.model.1"

#' Data: A small compendium of \emph{C. elegans} microarray data
#'
#' A \code{10,000 x 50} matrix containing Affymetrix microarray data gathered from
#' the GEO database for the model orgnaism C. elegans.  This matrix is a small subset
#' of the \code{GPL200.1386} matrix from the DEXDATA.Celegans
#' package intended for use in example scripts.  Rows of the matrix (probesets)
#' were chosen at random, but specific columns were selected in order to capture all of
#' the arrays from experiments with the fewest number of arrays, up to 50 arrays total.
#'
#' @format A numeric matrix with probesets in rows and arrays in columns.
#'
"GPL200.50.mini"

#' Data: A binary matrix associating probesets with Gene Ontology terms
#'
#' This is a binary probeset annotation matrix containing Affymetrix probeset IDs in rows
#' and Gene Ontology (GO) IDs in columns. A value of 1 in the matrix indicates that the
#' gene to which the probeset has been mapped has been annotated with the corresponding
#' GO term.  This matrix is a subset of the \code{GPL200.GO}
#' matrix from the DEXDATA.Celegans package and is intended for use in example scripts.
#' The matrix was constructed by first limiting the rows (probesets) to those found in
#' \code{\link[DEXICA]{GPL200.50.mini}}, then by selecting 1000 columns (GO IDs) at random
#' from those that contained at least 5 annoted probesets.
#'
#' @format A binary matrix with probeset IDs in rows and GO IDs in columns.
#'
"GPL200.GO.mini"
