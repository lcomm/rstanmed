# Global variables from non-standard evaluation
utils::globalVariables(c("am_intx", "bias", "ci_cov", "ci_width", 
                         "estimate", "estimator", "full_width", "n", 
                         "nder_truth", "transport"))

#' Process results from all registries and put in long format data frames
#' 
#' @param res_dir Result directory containing files named things like result_f5000_1p5.rds
#' @param ns Sample sizes of all registries in the directory
#' @param ttypes Transportability type suffixes. Defaults to 1p5 and 0_1p5
#' @return Named list with all results, all estimates, all biases,
#' all coverages, and all CI widths
#' @examples
#' \dontrun{
#' res <- compile_results(res_dir = "~/Dropbox/mediation-BSA/", 
#'                        ns = c(5000, 10000), ttypes = c("1p5", "0_1p5"))
#' }
#' @export
compile_results <- function(res_dir, ns = c(5000, 10000), ttypes = c("1p5", "0_1p5")) {
  i <- 0
  res <- list()
  long_width <- long_cov <- long_bias <- long_est <- list()
  suffixes <- c("_uc", "_dg", "_ix", "_gc", "_gf")
    
  for (rtype in c("f", "b")) {
    for (ttype in ttypes) {
      for (n in ns) {
        i <- i + 1
        fname <- paste0(res_dir, "result_", rtype, n, "_", ttype, ".rds")
        if (!file.exists(fname)) {
            message(paste0("Could not find result file ", fname))
            next
        }
        res[[i]] <- readRDS(file = fname)
        res[[i]]$n <- n
        res[[i]]$transport <- factor((ttype == "1p5") * 1, 
                                     levels = 0:1, 
                                     labels = c("No transportability",
                                                "Transportability"))
        res[[i]]$am_intx <- factor(res[[i]]$am_intx, 
                                   levels = 0:1, 
                                   labels = c("No interaction",
                                              "Exposure-mediator interaction"))
        id_vars <- c("seed", "u_ei", "am_intx", "transport", "nder_truth", "n")
        estimate_vars <- paste0("nder", suffixes)
        bias_vars <- paste0("bias", suffixes)
        cov_vars <- paste0("ci_cov", suffixes)
        width_vars <- paste0("ci_width", suffixes)
        long_est[[i]] <- res[[i]][, colnames(res[[i]]) %in% 
                                      c(id_vars, estimate_vars)]
        long_est[[i]] <- reshape2::melt(long_est[[i]], id.vars = id_vars,
                                        variable.name = "estimator",
                                        value.name = "estimate")
        long_bias[[i]] <- res[[i]][, colnames(res[[i]]) %in% 
                                       c(id_vars, bias_vars)]
        long_bias[[i]] <- reshape2::melt(long_bias[[i]], id.vars = id_vars,
                                         variable.name = "estimator",
                                         value.name = "bias")
        long_cov[[i]] <- res[[i]][, colnames(res[[i]]) %in% 
                                      c(id_vars, cov_vars)]
        long_cov[[i]] <- reshape2::melt(long_cov[[i]], id.vars = id_vars,
                                        variable.name = "estimator",
                                        value.name = "ci_cov")
        long_width[[i]] <- res[[i]][, colnames(res[[i]]) %in% 
                                        c(id_vars, width_vars)]
        long_width[[i]] <- reshape2::melt(long_width[[i]], id.vars = id_vars,
                                          variable.name = "estimator",
                                          value.name = "ci_width")
      }
    }
  }
    
  # Make one big long data set for each quantity
  full_est   <- as.data.frame(data.table::rbindlist(long_est))
  full_bias  <- as.data.frame(data.table::rbindlist(long_bias))
  full_cov   <- as.data.frame(data.table::rbindlist(long_cov))
  full_width <- as.data.frame(data.table::rbindlist(long_width))

  # Label levels for plots
  labels <- c("Naive", "DG", "IX", "BDF-SIM", "BDF-CF")
  full_est$estimator <- factor(full_est$estimator,
                               levels = paste0("nder", suffixes),
                               labels = labels)
  full_bias$estimator <- factor(full_bias$estimator,
                                levels = paste0("bias", suffixes),
                                labels = labels)
  full_cov$estimator <- factor(full_cov$estimator,
                               levels = paste0("ci_cov", suffixes),
                               labels = labels)
  full_width$estimator <- factor(full_width$estimator,
                                 levels = paste0("ci_width", suffixes),
                                 labels = labels)
  
  return(list(res = res, full_est = full_est, full_bias = full_bias,
              full_cov = full_cov, full_width = full_width))
  
}



#' Make the boxplot showing simulation estimates
#' 
#' Restricts to one level of u_ei
#' 
#' @param full_est Long data frame with estimates from all simulations. See 
#' \code{\link{compile_results}}.
#' @param u_ei Whether confounding is exposure-induced
#' @return ggplot2 object
#' @examples
#' \dontrun{
#' res <- compile_results(res_dir = "~/Dropbox/mediation-BSA/", 
#'                        ns = c(5000, 10000), ttypes = c("1p5", "0_1p5"))
#' png(filename = "~/Dropbox/mediation-BSA/sim_boxplot.png",
#'     width = 6, height = 6, units = "in", res = 300)
#' make_simulation_boxplot(full_est = res$full_est, u_ei = 1)
#' dev.off()
#' }
#' @export
make_simulation_boxplot <- function(full_est, u_ei) {
    
    # Subset
    est <- full_est[full_est$u_ei == u_ei, ]
    
    # Make plot
    p <- ggplot(data = est, aes(x = as.factor(n), y = estimate, fill = estimator)) +
      geom_boxplot(position=position_dodge(0.9), outlier.size = 0.7) +
      geom_hline(aes(yintercept = nder_truth, color = as.factor(nder_truth),
                     linetype = ""), color = grey.colors(1), show.legend = TRUE) +
      facet_grid(vars(am_intx), vars(transport)) +
      scale_linetype_manual(values = c("dotted")) +
      scale_fill_grey(start = 0.4, end = 0.95) +
      theme_bw() +
      lims(y = c(0.15, 0.42)) +
      labs(linetype = "Truth", 
           x = "Main sample size",
           y = "Estimated randomized natural direct effect",
           main = "Estimates from 100 data pair replicates",
           fill = "Estimator") +
      theme(legend.position = "bottom") +
      guides(fill = guide_legend(order = 2),
             linetype = guide_legend(order = 1))
    return(p)
}



#' Translate the main sample size $n$ into a string with the sample size pairs
#' 
#' Used for tables
#' 
#' @param n Main study sample size
#' @return String with (n, n / 10)
#' @export
prettify_main_n <- function(n){
  return(paste0("(", floor(n / 10), ", ", n,")"))
}



#' Make a table of all estimators bias within a given level of u_ei
#' 
#' @param full_bias Data frame containing long-form bias from replicates
#' @param u_ei Level of exposure-induced confounding for which to make the table
#' @param prettify Whether to prettify n printing to include both sample sizes
#' @return Data frame to be kabled
#' @examples 
#' \dontrun{
#' res <- compile_results(res_dir = "~/Dropbox/mediation-BSA/", 
#'                        ns = c(5000, 10000), ttypes = c("1p5", "0_1p5"))
#' knitr::kable(make_bias_table(res$full_bias, u_ei = 1, prettify = TRUE), 
#'              format = "latex", booktabs = TRUE, digits = 3, 
#'              caption = "Mean bias for the randomized natural direct effect in 
#'              naive, delta-gamma (DG) and interaction (IX) frequentist corrections, 
#'              simulation-based (BDF-SIM) and closed-form (BDF-CF) Bayesian data 
#'              fusion estimators, calculated in 100 replicates with exposure-induced 
#'              mediator-outcome confounding")
#' }
#' @export
make_bias_table <- function(full_bias, u_ei, prettify = FALSE) {
  a <- dplyr::group_by(.data = full_width,
                       estimator, n, u_ei, am_intx, transport)
  a <- dplyr::summarise(.data = a, bias = mean(bias)) 
  a <- a[a$u_ei == u_ei, -which(colnames(a) == "u_ei")]
  a$transport <- factor(a$transport, 
                        levels = c("No transportability","Transportability"),
                        labels = c("No", "Yes"))
  a$am_intx <- factor(a$am_intx, 
                      levels = c("No interaction", "Exposure-mediator interaction"),
                      labels = c("No", "Yes"))
  if (prettify) {
    a$ss <- prettify_main_n(a$n)
    b <- suppressMessages(reshape2::dcast(data = a, 
                                          transport + am_intx + ss ~ estimator,
                                          value.var = "bias"))
    colnames(b) <- c("Transportability", "Interaction", "Sample sizes",
                     c("Naive", "DG", "IX", "BDF-SIM", "BDF-CF"))
  } else {
    b <- suppressMessages(reshape2::dcast(data = a, 
                                          transport + am_intx + n ~ estimator,
                                          value.var = "bias"))
    colnames(b) <- c("Transportability", "Interaction", "Main study $n$",
                     c("Naive", "DG", "IX", "BDF-SIM", "BDF-CF"))
  }
  return(b)
}



#' Make a table of all estimators CI coverage within a given level of u_ei
#' 
#' @param full_cov Data frame containing long-form coverage from replicates
#' @param u_ei Level of exposure-induced confounding for which to make the table
#' @param prettify Whether to prettify n printing to include both sample sizes
#' @return Data frame to be kabled
#' @examples
#' \dontrun{
#' res <- compile_results(res_dir = "~/Dropbox/mediation-BSA/", 
#'                        ns = c(5000, 10000), ttypes = c("1p5", "0_1p5"))
#' knitr::kable(make_coverage_table(res$full_cov, u_ei = 1, prettify = TRUE),
#'                                  format = "latex", 
#'                                  booktabs = TRUE, digits = 3, 
#'                                  caption = "Coverage percentages for 95% confidence 
#'                                  and credible intervals for naive, delta-gamma 
#'                                  (DG) and interaction (IX) frequentist corrections, 
#'                                  simulation-based (BDF-SIM) and closed-form (BDF-CF) 
#'                                  Bayesian data fusion estimators, calculated in 100
#'                                  replicates with exposure-induced mediator-outcome 
#'                                  confounding")
#' }
#' @export
make_coverage_table <- function(full_cov, u_ei, prettify = FALSE) {
  a <- dplyr::group_by(.data = full_width,
                       estimator, n, u_ei, am_intx, transport)
  a <- dplyr::summarise(.data = a, ci_cov = mean(ci_cov) * 100) 
  a <- a[a$u_ei == u_ei, -which(colnames(a) == "u_ei")]
  a$transport <- factor(a$transport, 
                        levels = c("No transportability","Transportability"),
                        labels = c("No", "Yes"))
  a$am_intx <- factor(a$am_intx, 
                      levels = c("No interaction", "Exposure-mediator interaction"),
                      labels = c("No", "Yes"))
  if (prettify) {
    a$ss <- prettify_main_n(a$n)
    b <- suppressMessages(reshape2::dcast(data = a, 
                        transport + am_intx + ss ~ estimator,
                        value.var = "ci_cov"))
    colnames(b) <- c("Transportability", "Interaction", "Sample sizes",
             c("Naive", "DG", "IX", "BDF-SIM", "BDF-CF"))
  } else {
    b <- suppressMessages(reshape2::dcast(data = a, 
                        transport + am_intx + n ~ estimator,
                        value.var = "ci_cov"))
    colnames(b) <- c("Transportability", "Interaction", "Main study $n$",
             c("Naive", "DG", "IX", "BDF-SIM", "BDF-CF"))
  }
  return(b)
}



#' Make a table of all estimators CI coverage within a given level of u_ei
#' 
#' @param full_width Data frame containing long-form CI widths from replicates
#' @param u_ei Level of exposure-induced confounding for which to make the table
#' @param prettify Whether to prettify n printing to include both sample sizes
#' @return Data frame to be kabled
#' @examples 
#' \dontrun{
#' res <- compile_results(res_dir = "~/Dropbox/mediation-BSA/", 
#'                        ns = c(5000, 10000), ttypes = c("1p5", "0_1p5"))
#' knitr::kable(make_width_table(res$full_width, u_ei = 1, prettify = TRUE),
#'              format = "latex", 
#'              booktabs = TRUE, digits = 3, 
#'              caption = "Widths of 95% confidence and credible intervals for
#'              naive, delta-gamma (DG) and interaction (IX) frequentist corrections, 
#'              simulation-based (BDF-SIM) and closed-form (BDF-CF) Bayesian data 
#'              fusion estimators, calculated in 100 replicates with exposure-induced 
#'              mediator-outcome confounding")
#' }
#' @export
make_width_table <- function(full_width, u_ei, prettify = FALSE) {
  a <- dplyr::group_by(.data = full_width,
                       estimator, n, u_ei, am_intx, transport)
  a <- dplyr::summarise(.data = a, ci_width = mean(ci_width))
  a <- a[a$u_ei == u_ei, -which(colnames(a) == "u_ei")]
  a$transport <- factor(a$transport, 
                        levels = c("No transportability","Transportability"),
                        labels = c("No", "Yes"))
  a$am_intx <- factor(a$am_intx, 
                      levels = c("No interaction", "Exposure-mediator interaction"),
                      labels = c("No", "Yes"))
  if (prettify) {
    a$ss <- prettify_main_n(a$n)
    b <- suppressMessages(reshape2::dcast(data = a, 
                                          transport + am_intx + ss ~ estimator,
                                          value.var = "ci_width"))
    colnames(b) <- c("Transportability", "Interaction", "Sample sizes",
                     c("Naive", "DG", "IX", "BDF-SIM", "BDF-CF"))
  } else {
    b <- suppressMessages(reshape2::dcast(data = a, 
                                          transport + am_intx + n ~ estimator,
                                          value.var = "ci_width"))
    colnames(b) <- c("Transportability", "Interaction", "Main study $n$",
                     c("Naive", "DG", "IX", "BDF-SIM", "BDF-CF"))
  }
  return(b)
}


