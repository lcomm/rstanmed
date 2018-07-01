#' Break up existing chunks of registry ids
#' 
#' @param reg batchtools registry object
#' @param n.chunks Number of chunks to be submitted
#' @param chunk.size Size of chunks to be submitted
#' @return Chunked ids to be passed to batchtools::submitJobs
#' @export
chunk_registry <- function(reg, n.chunks = NULL, chunk.size = NULL) {
  ids <- batchtools::getJobTable(reg = reg)
  ids$chunk = batchtools::chunk(ids$job.id, 
                                n.chunks = n.chunks, 
                                chunk.size = chunk.size)
  return(ids)
}


#' Function to submit many BDF simulation study replicates using a batchtools
#' registry
#' 
#' @param registry batchtools registry object
#' @param transport TRUE/FALSE for whether to make effects transportable
#' @param R Number of replicates
#' @param clear_existing Whether to clear existing registry first
#' @param u_ei Exposure-induced flag values
#' @param am_intx Exposure-mediator flag values
#' @param n Values for (big) sample sizes
#' @param result_type Output type, from \code{\link{run_bdf_replicate}}
#' @param iter Number of MCMC iterations
#' @param chains Number of MCMC chains
#' @param chunk.size How many jobs should be chunked together
#' @param time_each Number of minutes for each job at the smallest n
#' @param memory Memory to allocate at the smallest n
#' @param max.concurrent.jobs Maximum number of jobs at the same time
#' @return None; jobs will be submitted and update in registry
#' @export
submit_bdf_jobs <- function(registry, transport, R, 
                            clear_existing = FALSE,
                            u_ei = 0:1, am_intx = 0:1, 
                            n = c(5000, 10000),
                            result_type = "processed", 
                            iter = 2000, chains = 4,
                            chunk.size = 4,
                            time_each = 120,
                            memory = 4000,
                            max.concurrent.jobs = 2000) {
  if (clear_existing) {
    batchtools::clearRegistry(reg = registry)
  } else {
    if (NROW(batchtools::getJobTable(reg = registry)) > 0) {
      stop("Registry has existing jobs but clear_existing is FALSE")  
    }
  }
  
  prog_opt <- getOption("batchtools.progress")
  options(batchtools.progress = FALSE)
  
  if (transport) {
    strengths <- rep(1.5, 4)
  } else {
    strengths <- c(1.5, 1.5, 0, 0)
  }
  
  # Make job
  args <- data.table::CJ(u_ei = u_ei, am_intx = am_intx, seed = 1:R,
                         yu_strength = strengths[1], 
                         mu_strength = strengths[2], 
                         small_yu_strength = strengths[3], 
                         small_mu_strength = strengths[4],
                         n = n)
  batchtools::batchMap(fun = run_bdf_replicate, 
         args = args, 
         more.args = list(prior_type = "dd", 
                          params = NULL,
                          small_params = NULL,
                          result_type = result_type,
                          iter = iter, chains = chains), 
         reg = registry)
  
  # Decide chunking, keeping n homogeneous within chunks
  # ids <- cbind(getJobTable(reg = registry)[, "job.id"], n = args$n)
  # ids[, n := args$n]
  # ids[, chunk := chunk(job.id, chunk.size = chunk.size), by = "n"]
  # ids[, chunk := .GRP, by = c("n", "chunk")]
  
  # Submit job chunks grouped by sample size n
  # for (n_val in unique(n)) {
    n_ratio <- 1 #round(n_val / min(n))
    walltime <- 60 * n_ratio * time_each * chunk.size
    memory <- n_ratio * memory
    batchtools::submitJobs(ids = chunk_registry(reg = registry,
                                                chunk.size = chunk.size),
                           # ids = ids[ids$n == n_val, -"n"],
                           reg = registry,
                           resources = list(walltime = walltime,
                                            memory = memory,
                                            max.concurrent.jobs = 
                                            max.concurrent.jobs))
  # }
  
  # Reset option
  options(batchtools.progress = prog_opt)
}
