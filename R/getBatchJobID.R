getBatchJobID = function(job.index.ID = c("SLURM_ARRAY_TASK_ID", "SGE_TASK_ID",
                                          "PBS_ARRAYID", "LSB_JOBINDEX"), default = 1) {
  for(ID in job.index.ID) {
    job.ID = as.numeric(Sys.getenv(ID))
    if(!is.na(job.ID)) return(job.ID)
  }

  # No valid job ID found.
  warning("Could not determine batch job ID, returning default value: ", default)
  return(default)
}
