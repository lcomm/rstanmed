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



#' Conditionally clear a registry if requested
#' 
#' @param clear_existing Whether to clear registry
#' @param registry Batchtools registry to clear (or not)
#' @return None (if registry empty), or error message (if not empty and not
#' requested to clear). Used for side effects.
#' @export
process_clear_registry <- function(clear_existing, registry) {
  if (clear_existing) {
    batchtools::clearRegistry(reg = registry)
  } else {
    if (NROW(batchtools::getJobTable(reg = registry)) > 0) {
      stop("Registry has existing jobs but clear_existing is FALSE")  
    }
  }
}



#' Function to submit many BDF simulation study replicates using a batchtools
#' registry
#' 
#' @param registry batchtools registry object
#' @param transport TRUE/FALSE for whether to make effects transportable
#' @param seed Seed vector (\code{seed = 1:R} does R replicates)
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
submit_bdf_jobs <- function(registry, transport, seed, 
                            clear_existing = FALSE,
                            u_ei = 0:1, am_intx = 0:1, 
                            n = c(5000, 10000),
                            result_type = "processed",
                            iter = 2000, chains = 4,
                            chunk.size = 4,
                            time_each = 120,
                            memory = 4000,
                            n_ratio = 10,
                            max.concurrent.jobs = 2000) {
  
  process_clear_registry(clear_existing, registry)
  
  prog_opt <- getOption("batchtools.progress")
  options(batchtools.progress = FALSE)
  
  # Make job
  strengths <- convert_transport_to_strengths(transport)
  args <- data.table::CJ(u_ei = u_ei, am_intx = am_intx, seed = seed,
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
                          n_ratio = n_ratio,
                          iter = iter, chains = chains), 
         reg = registry)
  
  walltime <- 60 * time_each * chunk.size
  batchtools::submitJobs(ids = chunk_registry(reg = registry,
                                              chunk.size = chunk.size),
                         reg = registry,
                         resources = list(walltime = walltime,
                                          memory = memory,
                                          max.concurrent.jobs = 
                                            max.concurrent.jobs))
  
  # Reset option
  options(batchtools.progress = prog_opt)
}



#' Function to submit many BDF simulation study replicates using a batchtools
#' registry
#' 
#' @param registry batchtools registry object
#' @param transport TRUE/FALSE for whether to make effects transportable
#' @param seed Seed vector (\code{seed = 1:R} does R replicates)
#' @param clear_existing Whether to clear existing registry first
#' @param u_ei Exposure-induced flag values
#' @param am_intx Exposure-mediator flag values
#' @param n Values for (big) sample sizes
#' @param result_type Output type, from \code{\link{run_frequentist_replicate}}
#' @param n_ratio 
#' @param B Number of bootstrap samples for CI
#' @param chunk.size How many jobs should be chunked together
#' @param time_each Number of minutes for each job at the smallest n
#' @param memory Memory to allocate at the smallest n
#' @param max.concurrent.jobs Maximum number of jobs at the same time
#' @return None; jobs will be submitted and update in registry
#' @export
submit_frequentist_jobs <- function(registry, transport, seed, 
                                    clear_existing = FALSE,
                                    u_ei = 0:1, am_intx = 0:1, 
                                    n = c(5000, 10000),
                                    result_type = "processed", 
                                    n_ratio = 10,
                                    B = 1000,
                                    chunk.size = 20,
                                    time_each = 15,
                                    memory = 4000,
                                    max.concurrent.jobs = 2000) {
  
  process_clear_registry(clear_existing, registry)
  
  prog_opt <- getOption("batchtools.progress")
  options(batchtools.progress = FALSE)
  
  # Make job
  strengths <- convert_transport_to_strengths(transport)
  args <- data.table::CJ(u_ei = u_ei, am_intx = am_intx, seed = seed,
                         yu_strength = strengths[1], 
                         mu_strength = strengths[2], 
                         small_yu_strength = strengths[3], 
                         small_mu_strength = strengths[4],
                         n = n)
  batchtools::batchMap(fun = run_frequentist_replicate, 
                       args = args, 
                       more.args = list(params = NULL,
                                        small_params = NULL,
                                        result_type = result_type,
                                        n_ratio = n_ratio,
                                        B = B), 
                       reg = registry)
  
  walltime <- 60 * time_each * chunk.size
  batchtools::submitJobs(ids = chunk_registry(reg = registry,
                                              chunk.size = chunk.size),
                         reg = registry,
                         resources = list(walltime = walltime,
                                          memory = memory,
                                          max.concurrent.jobs = 
                                            max.concurrent.jobs))
  
  # Reset option
  options(batchtools.progress = prog_opt)
}



#' Make a shell data frame to contain job parameters from registry output
#' 
#' @param registry Batchtools registry containing jobs
#' @param registry_type Whether registry has frequentist ("freq") or Bayesian ("bdf")
#' results
#' @return Data frame containing job parameters and empty columns for results
#' @export
make_simregdf_shell <- function(registry, registry_type = c("freq", "bdf")) {
  job_parl <- getJobPars(reg = registry)$job.pars
  if (registry_type == "freq") {
    df <- data.frame(job_id = getJobPars(reg = registry)$job.id,
                     seed = sapply(job_parl, function(x) x$seed),
                     u_ei = sapply(job_parl, function(x) x$u_ei),
                     am_intx = sapply(job_parl, function(x) x$am_intx),
                     n = sapply(job_parl, function(x) x$n),
                     yu_strength = sapply(job_parl, function(x) x$yu_strength),
                     mu_strength = sapply(job_parl, function(x) x$mu_strength),
                     small_yu_strength = sapply(job_parl, 
                                                function(x) x$small_yu_strength),
                     small_mu_strength = sapply(job_parl, 
                                                function(x) x$small_mu_strength),
                     nder_truth = NA,
                     nder_uc = NA,
                     nder_dg = NA,
                     nder_ix = NA,
                     bias_uc = NA,
                     bias_dg = NA,
                     bias_ix = NA,
                     ci_cov_uc = NA,
                     ci_cov_dg = NA,
                     ci_cov_ix = NA,
                     ci_width_uc = NA,
                     ci_width_dg = NA,
                     ci_width_ix = NA)
  } else if (registry_type == "bdf") {
    df <- data.frame(job_id = getJobPars(reg = registry)$job.id,
                     seed = sapply(job_parl, function(x) x$seed),
                     u_ei = sapply(job_parl, function(x) x$u_ei),
                     am_intx = sapply(job_parl, function(x) x$am_intx),
                     n = sapply(job_parl, function(x) x$n),
                     nder_truth = NA,
                     nder_gc = NA,
                     nder_gf = NA,
                     bias_gc = NA,
                     bias_gf = NA,
                     ci_cov_gc = NA,
                     ci_cov_gf = NA,
                     ci_width_gc = NA,
                     ci_width_gf = NA)
  } else {
    stop("Unknown registry type!")
  }
  return(df)
}



#' Process registry of frequentist results
#' 
#' @param registry Registry containing jobs
#' @return Data frame containing parameter conditions and results
#' @export
combine_frequentist_reg_results <- function(registry) {
  shell <- make_simregdf_shell(registry = registry, registry_type = "freq")
  for (i in 1:NROW(shell)) {
    rdsname <- paste0(path.expand(registry$file.dir), "/results/", i, ".rds")
    is_done <- (!is.na(getJobStatus(id = i, reg = registry)$done)) && 
               (file.exists(rdsname))
    if (is_done) {
      job_res <- batchtools::loadResult(id = i, reg = registry)
      shell[i, "nder_truth"] <- job_res$truth_nder
      for (suffix in c("uc", "dg", "ix")) {
        estimate_name <- paste0("nder_", suffix)
        bias_name     <- paste0("bias_", suffix)
        cov_name      <- paste0("ci_cov_", suffix)
        width_name    <- paste0("ci_width_", suffix)
        shell[i, estimate_name] <- job_res$estimate[[suffix]]
        shell[i, bias_name]     <- job_res$bias[[suffix]]
        shell[i, cov_name]      <- job_res$ci_cov[[suffix]]
        shell[i, width_name]    <- job_res$ci_width[[suffix]]
      }
    }
  }
  return(shell)
}



#' Process registry of Bayesian Data Fusion results
#' 
#' @param registry Registry containing jobs
#' @return Data frame containing parameter conditions and results
#' @export
combine_bdf_reg_results <- function(registry) {
  shell <- make_simregdf_shell(registry = registry, registry_type = "bdf")
  for (i in 1:NROW(shell)) {
    rdsname <- paste0(path.expand(registry$file.dir), "/results/", i, ".rds")
    is_done <- (!is.na(getJobStatus(id = i, reg = registry)$done)) && 
      (file.exists(rdsname))
    if (is_done) {
      job_res <- batchtools::loadResult(id = i, reg = registry)
      shell[i, "nder_truth"] <- job_res$truth_nder
      for (suffix in c("gc", "gf")) {
        estimate_name <- paste0("nder_", suffix)
        bias_name     <- paste0("bias_", suffix)
        cov_name      <- paste0("ci_cov_", suffix)
        width_name    <- paste0("ci_width_", suffix)
        shell[i, estimate_name] <- job_res$estimate[[suffix]]
        shell[i, bias_name]     <- job_res$bias[[suffix]]
        shell[i, cov_name]      <- job_res$ci_cov[[suffix]]
        shell[i, width_name]    <- job_res$ci_width[[suffix]]
      }
    }
  }
  return(shell)
}
