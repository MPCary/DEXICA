################################
### Create class definitions ###
################################


# Required to allow w.init to be either a number or NULL:
setClassUnion("numericOrNULL", c("numeric", "NULL"))

#' An S4 class for defining DexBatch parameter ranges
#'
#'@description A ParameterSet object defines a range of parameter values to
#'test in a DexBatch (see \code{\link[DEXICA]{DexBatch-class}}).
#'
#'@usage ## Constructor:
#'parameterSet(n.comp = 1, center.cols = c(TRUE), scale.cols = c(FALSE),
#'center.rows = c(FALSE), scale.rows = c(FALSE),
#'alg.typ = "parallel", fun = "logcosh", alpha = 1,
#'maxit = 200, tol = 1e-04, w.init = NULL, max.attempts = 3, verbose = FALSE)
#'
#'@slot n.comp an integer or vector specifying the number of components to extract
#'@slot alg.typ,fun,alpha,maxit,tol,w.init,max.attempts,verbose See \code{\link[DEXICA]{dexFastICA}} documentation
#'@slot center.cols,scale.cols,center.rows,scale.rows,simultaneous See \code{\link[DEXICA]{preprocessMatrix}} documentation
#'
#'@details When creating a DexBatch object, a ParameterSet may
#'be supplied (in the \code{params} slot) in order to specify a range of parameter
#'values to test.  Each slot in ParameterSet may accept either a single value or a
#'range of values.  If multiple slots have multiple values, the ParameterSet will
#'encompass each possible combination of individual slot values.
#'
#'Any slot values not specified during object construction are given their default
#'values (see constructor), e.g., calling \code{parameterSet()}
#'results in a ParameterSet with a single set of parameter values.  The
#'\code{\link[DEXICA]{getParamSetCount}} method returns the total number of
#'parameter sets defined in a ParameterSet.
#'
#'@export
parameterSet = setClass(
  # An S4 object is nessary here in order to prevent batch jobs from being
  # called with incomplete or non-sense parameter lists.  A ParameterSet
  # object is meant to define the range of parameters that a batch of ICA
  # runs will explore.  Calling getParams(m) will return the specific set
  # of parameters for job m.

  "ParameterSet",

  # Define slots
  slots = c(
    n.comp = "numeric",
    center.cols = "logical",
    scale.cols = "logical",
    center.rows = "logical",
    scale.rows = "logical",
    alg.typ = "character",
    fun = "character",
    alpha = "numeric",
    maxit = "numeric",
    tol = "numeric",
    w.init = "numericOrNULL",
    max.attempts = "numeric",
    verbose = "logical",
    param.space = "data.frame"
  )
)

setMethod(f = "initialize", signature = "ParameterSet",
          function(.Object, n.comp = 1,
                   center.cols = c(TRUE),
                   scale.cols = c(FALSE),
                   center.rows = c(FALSE),
                   scale.rows = c(FALSE),
                   alg.typ = "parallel",
                   fun = "logcosh",
                   alpha = 1,
                   maxit = 200,
                   tol = 1e-04,
                   w.init = NULL,
                   max.attempts = 3,
                   verbose = FALSE, ...) {

  # Set parameters
  .Object@n.comp = n.comp
  .Object@center.cols = center.cols
  .Object@scale.cols = scale.cols
  .Object@center.rows = center.rows
  .Object@scale.rows = scale.rows
  .Object@alg.typ = alg.typ
  .Object@fun = fun
  .Object@alpha = alpha
  .Object@maxit = maxit
  .Object@tol = tol
  .Object@max.attempts = max.attempts
  .Object@verbose = verbose

  # Pick a random seed if w.init is NULL
  if(is.null(w.init)) w.init = round(runif(1, min = 0, max = 10^6))
  .Object@w.init = w.init

  # Build table of parameter sets to test
  param.space = expand.grid(n.comp = n.comp, center.cols = center.cols,
                          scale.cols = scale.cols, center.rows = center.rows,
                          scale.rows = scale.rows,
                          alg.typ = alg.typ, fun = fun,
                          alpha = alpha, maxit = maxit, tol = tol, w.init = w.init,
                          max.attempts = max.attempts, verbose = verbose,
                          stringsAsFactors = FALSE)
  .Object@param.space = param.space

  # Check validity
  validObject(.Object)

  # Return object
  .Object
})

validParameterSet = function(object) {
  if(!(is.numeric(object@w.init))) {
      return("w.init must be numeric.")
    }
  TRUE
}

setValidity("ParameterSet", method = validParameterSet)


################################
### Create generic functions ###
################################

#' Get the number of parameter sets defined in a ParameterSet object
#'
#' @description This function returns the number of individual parameter sets defined in a
#' ParameterSet object.
#'
#' @param object Object to get parameter set count for.
#'
#' @examples
#' p = parameterSet(n.comp = c(10,20,30))
#' getParamSetCount(p)
#'
setGeneric(name = "getParamSetCount", def = function(object) {
  standardGeneric("getParamSetCount")
})

#' Get a single parameter set defined in a ParameterSet object
#'
#' @description This function returns a single parameter sets defined in a
#' ParameterSet object.
#'
#' @param object Object to get parameter set from.
#' @param j An integer indicating which parameter set in object to return
#'
#' @examples
#' p = parameterSet(n.comp = c(10,20,30))
#' getOneParamSet(p, 1)
#'
setGeneric(name = "getOneParamSet", def = function(object, j) {
  standardGeneric("getOneParamSet")
})

#' Get all parameter sets defined in a ParameterSet object
#'
#' @description This function returns all parameter sets defined in a
#' ParameterSet object.
#'
#' @param object Object to get parameter sets from.
#'
#' @examples
#' p = parameterSet(n.comp = c(10,20,30))
#' getAllParamSets(p)
#'
setGeneric(name = "getAllParamSets", def = function(object) {
  standardGeneric("getAllParamSets")
})

#' Reset the random seeds in a ParameterSet object
#'
#' @description This function resets the random seeds defined in a ParameterSet
#' (specifically, in the \code{w.init} slot) to random values if new.seeds is NULL,
#' or to the values specified in new.seeds otherwise.
#'
#' @param object Object in which to reset random seeds
#' @param new.seeds Numeric vector indicating new random seeds
#'
#' @details If new.seeds is not NULL, it must be a numeric vector with a length
#' equal to the number of parameter sets defined in \code{object} (this value
#' can be determined with the \code{\link[DEXICA]{getParamSetCount}} method.)
#'
#' @examples
#' p = parameterSet(n.comp = c(10,20,30))
#' getAllParamSets(p)
#'
setGeneric(name = "reseedParamSets", def = function(object, new.seeds = NULL) {
  standardGeneric("reseedParamSets")
})


######################
### Create methods ###
######################

#' @describeIn ParameterSet Get the number of parameter sets defined in a ParameterSet object (see
#' \code{\link[DEXICA]{getParamSetCount}})
setMethod(f = "getParamSetCount", signature = "ParameterSet",
          function(object) {
            nrow(object@param.space)
          })

#' @describeIn ParameterSet Get a single parameter set defined in a ParameterSet object (see
#' \code{\link[DEXICA]{getOneParamSet}})
setMethod(f = "getOneParamSet", signature = "ParameterSet",
          function(object, j) {
            n = getParamSetCount(object)
            if((j >= 1) & (j <= n)) {
              return(object@param.space[j, ])
              } else {
              stop("j must be a number from 1 to ", n)
            }
})

#' @describeIn ParameterSet Get all parameter sets defined in a ParameterSet object (see
#' \code{\link[DEXICA]{getAllParamSets}})
setMethod(f = "getAllParamSets", signature = "ParameterSet",
          function(object) {
            object@param.space
          })

#' @describeIn ParameterSet Reset the random seeds in a ParameterSet object (see
#' \code{\link[DEXICA]{reseedParamSets}})
setMethod(f = "reseedParamSets", signature = "ParameterSet",
          function(object, new.seeds = NULL) {
            if(is.null(new.seeds)) {
              n = getParamSetCount(object)
              new.seeds = round(runif(n, min = 0, max = 10^6))
            } else {
              # Set to user-supplied values
              if(!is.numeric(new.seeds)) stop("new.seeds must be numeric")
              if(!(length(new.seeds) == getParamSetCount(object))) {
                stop("length of new.seeds must match number of parameter sets")
              }
            }
            object@param.space[, 'w.init'] = new.seeds
            object
          })
