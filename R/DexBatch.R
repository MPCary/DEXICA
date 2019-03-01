#' @include ParameterSet.R

################################
### Create class definitions ###
################################

# To do: make a method for estimating the time per job and total run time

# Required to allow package and file.path to be either character or NULL:
setClassUnion("ListOrMatrix", c("list", "matrix"))
setClassUnion("characterOrNULL", c("character", "NULL"))

#' An S4 class for testing multiple compendia and/or parameters
#'
#'@description The DexBatch class helps users test the effect of different
#'parameter options and/or input compendia on gene module prediction.
#'
#'@usage ## Constructor:
#'dexBatch(compen, annmats, params = NULL, output = NULL)
#'
#'@slot compen A matrix or list of matrices containing the compendia to be tested
#'@slot annmats A matrix or list of matrices containing the annotation matrices to be used
#'in module evaltation (see \code{\link[DEXICA]{checkModules}})
#'@slot params A ParameterSet object containing desired parameter ranges to be
#'tested (see \code{\link[DEXICA]{parameterSet}})
#'@slot output Character string name of output file
#'
#'@details A DexBatch object defines a set of jobs.  Each job first predicts a set of
#'gene modules, then evaluates their quality.  Jobs within a DexBatch differ from
#'each other in some way; the purpose of a DexBatch is to test the effect of
#'these differences on the quality of the resulting modules.
#'
#'Individual jobs in a DexBatch are executed by calling \code{\link[DEXICA]{runJob}}.
#'While each job could be run sequentially on a single machine, it is usually
#'more practical when testing a large number of parameters to execude the jobs in
#'parallel (see Use Cases in the DEXICA vignette)
#'
#'If params is NULL, a single parameter set will be supplied using default values for
#'the \code{\link[DEXICA]{parameterSet}} constructor.  If output is NULL, results will
#'be directed to \code{stdout()}.
#'
#'@export
dexBatch = setClass("DexBatch",
                      slots = c("compen" = "list",
                                "annmats" = "list",
                                "params" = "ParameterSet",
                                "output" = "characterOrNULL",
                                ".job.table" = "data.frame"
                                )
)

setMethod(f= "initialize", signature = "DexBatch",
          function(.Object, compen = NULL, annmats = NULL, params = NULL, output = NULL) {
            .Object@compen = compen
            .Object@annmats = annmats
            if(is.null(params)) {
              .Object@params = parameterSet()
            } else {
              .Object@params = params
            }

            # Create connection (if needed)
            if(is.null(output)) {
              .Object@output = "stdout"
            } else {
              .Object@output = output
            }

            # Create job table
            .Object@.job.table = expand.grid(names(.Object@compen),
                                             c(1:getParamSetCount(.Object@params)))
            colnames(.Object@.job.table) = c("compendium", "parameter_set")

            # Check validity
            validObject(.Object)

            # Return .Object
            .Object
          })

validDexBatch = function(object) {
  if(is.null(object@compen)) return("one or more compendia must be specified.")
  if(is.null(object@annmats)) return("one or more annotation matrices must be specified.")
  if(is.null(names(object@compen))) return("compendia must be provided in a named list.")
  if(is.null(names(object@annmats))) return("annotation matrices must be provided in a named list.")

  # Check names length and uniqueness
  if(any(names(object@compen) == "")) return("all compendia must be named")
  if(any(names(object@annmats) == "")) return("all annotation matrices must be named")
  if(length(unique(names(object@compen))) != length(object@compen)) {
    return("all compendia must have unique names")
  }
  if(length(unique(names(object@annmats))) != length(object@annmats)) {
    return("all annotation matrices must have unique names")
  }

  # Test output connection
  if(object@output != "stdout") {
    con = file(object@output)
    open(con, open = "a")
    if(isOpen(con)) {
      # Connection opened successfully
      close(con)
    } else {
      return("could not open output file connection")
    }
  }
  TRUE
}

setValidity("DexBatch", method = validDexBatch)


################################
### Create generic functions ###
################################


#' Count the number of jobs defined in a DexBatch object
#'
#' @description This function returns the expected number of individual jobs defined in a
#' DexBatch object.  The formula it uses is \code{c * p}, where c is the number of
#' compendia and p is the number of different parameter sets to be tested.
#'
#' @param object Object to get job count for.
#'
#' @examples
#' x = matrix(rnorm(10 * 10), 10, 10)
#' y = matrix(rnorm(10 * 10), 10, 10)
#' a1 = matrix(sample(c(0,1), 100, replace = TRUE))
#' db = dexBatch(list(x,y), list(a1))
#' countJobs(db)
#'
setGeneric("countJobs", def = function(object) {
  standardGeneric("countJobs")
})

#' Run a single job defined in a DexBatch object
#'
#' @description This function executes a single job defined in a DexBatch object and writes
#' the results to the output connection defined therein.
#'
#' @param object DexBatch object that defines job parameters.
#' @param j Job number to run
#'
#' @examples
#' x = matrix(rnorm(10 * 10), 10, 10)
#' y = matrix(rnorm(10 * 10), 10, 10)
#' a1 = matrix(sample(c(0,1), 100, replace = TRUE))
#' db = dexBatch(list(x,y), list(a1))
#' runJob(db, 1)
#'
setGeneric("runJob", def = function(object, j) {
  standardGeneric("runJob")
})

######################
### Create methods ###
######################

#' @describeIn DexBatch Count the number of jobs defined in a DexBatch object (see
#' \code{\link[DEXICA]{countJobs}})
setMethod(f = "countJobs", signature = "DexBatch",
          function(object) {
            return(nrow(object@.job.table))
          })

#' @describeIn DexBatch Run a single job defined in a DexBatch object (see
#' \code{\link[DEXICA]{runJob}})
setMethod(f = "runJob", signature = "DexBatch",
          function(object, j) {

            # Get compendium and parameter set numbers
            my.compen.num = object@.job.table[j, "compendium"]
            my.paramSet.num = object@.job.table[j, "parameter_set"]

            # Get compendium and parameter set
            x = object@compen[[my.compen.num]]
            name.of.x = names(object@compen)[[my.compen.num]]
            my.params = getOneParamSet(object@params, my.paramSet.num)

            # Check verbosity
            verbose = my.params$verbose

            # Generate description of this run for results
            results = c("job.num" = j, "compen" = name.of.x, my.params)

            # Preprocess compendium
            if(verbose) message("Preprocessing compendium...")
            x.pp = preprocessMatrix(x, center.cols = my.params$center.cols,
                                    scale.cols = my.params$scale.cols,
                                    center.rows = my.params$center.rows,
                                    scale.rows = my.params$scale.rows,
                                    verbose = my.params$verbose)

            # Predict modules
            if(verbose) message("Predicting modules...")
            m = predictModules(x.pp, n.comp = my.params$n.comp,
                               alg.typ = my.params$alg.typ,
                               fun = my.params$fun,
                               alpha = my.params$alpha,
                               maxit = my.params$maxit,
                               tol = my.params$tol,
                               w.init = my.params$w.init,
                               max.attempts = my.params$max.attempts,
                               verbose = my.params$verbose)

            # Add module prediction output to results
            results = c(results, "ica.attempts" = m$attempts,
                        "ica.iterations" = m$iterations,
                        "ica.actual.tol" = m$actual.tol,
                        "ica.min.obs.tol" = m$min.obs.tol,
                        "ica.converged" = m$converged)

            # Get module definitions
            if(verbose) message("Partitioning modules into sets...")
            S.vp = partition(m$S, method = "ann") # Variable partitioning
            S.vp.split = splitMatrix(S.vp)

            # Check modules
            for(a in names(object@annmats)) {
              if(verbose) message("Checking modules for ", a, " enrichment...")
              # Run check
              my.check = checkModules(S.vp.split, as.matrix(object@annmats[[a]]))

              # Generate results for this annotation matrix
              names(my.check) = paste(a, names(my.check), sep = ".")
              results = c(results, my.check)
              if(verbose) message("done.")
            }

            # Write output
            if(verbose) message("Writing results to output connection...")
            if(object@output == "stdout") {
              write(unlist(results), file = stdout(), ncolumns = length(results),
                    append = TRUE, sep = "\t")
            } else {
              # Writing to a file
              # The user may have deleted the output file (e.g., from a previous)
              # run, so create it if needed.
              con = file(object@output)
              open(con, open = "a")
              if(isOpen(con)) {
                # Connection opened successfully
              } else {
                stop("could not open output file connection")
              }

              # If the file is empty, first write the column names
              if (file.info(object@output)$size == 0) {
                write(names(results), file= con, ncolumns = length(results),
                      append = TRUE, sep = "\t")
              }
              # Write run results
              write(unlist(results), file= con, ncolumns = length(results),
                         append = TRUE, sep = "\t")
              close(con)
            } # End write to file
          } # End runJob function
)



























# dexBatch = function(job_num, compen, annot, ParameterSet, results.file = "dexBatch.results.txt", p = 0.05) {
#   # Get parameters for this job and load the x matrix
#   params = getJobParameters(ParameterSet, j)
#   results = c(params, "p.value.threshold" = p) # Named list
#   x = data(eval(params$x))
#
#   # Preprocess matrix
#   x.pp = preprocessMatrix(x, center.cols = params$center.cols, scale.cols = params$scale.cols,
#                           center.rows = params$center.rows, scale.rows = params$scale.rows,
#                           verbose = params$verbose)
#
#   # Generate a random seed and record it
#   s = runif(1)
#   results = c(results, "random.seed" = s)
#   set.seed(s)
#
#   # Predict modules
#   m = predictModules(x.pp, n.comp = params$n.comp, alg.typ = params$alg.typ, fun = params$fun, alpha = params$alpha,
#                      maxit = params$maxit, tol = params$tol, w.init = params$w.init,
#                      max.attempts = params$max.attempts, verbose = params$verbose)
#
#   # Check modules
#   for(a in AnnotationMatrixList) {
#     # Run check
#     my.check = checkModules(m, local(get(load(a))))
#
#     # Generate named results
#     names(my.check) = paste(names(a), names(my.check), sep = ".")
#     results = c(results, my.check)
#   }
#
#   # Write results to file
#   write(results, file = file, append = TRUE, sep = "\t")
#
# }
