#' Simple binary (M) binary (Y) mediation analysis in Stan
#'
#' @param outcomes Named list of outcomes m, u, and y
#' @param designs Named list of design matrices for m, u, and y
#' @param ... Additional parameters to pass to `rstan::sampling`
#' @return an object of class `stanfit` returned by `rstan::sampling`
#' @export
bin_bin_stan <- function(outcomes, designs, ...) {
  out <- rstan::sampling(stanmodels$bin_bin, 
                         data=list(N = length(outcomes$y), 
                                   P_m = ncol(designs$x_m) + 1,
                                   P_y = ncol(designs$x_y) + 1,
                                   P_y = ncol(designs$x_u),
                                   y = outcomes$y,
                                   m = outcomes$m,
                                   u = outcomes$u,
                                   x_m = designs$x_m,
                                   x_y = designs$x_y,
                                   x_u = designs$x_u), ...)
  return(out)
  
}
