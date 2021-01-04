#' @title
#' Check error
#'
#' @param x return from try() function
#'
#' @return boolean indicating if the call resulted in error
#'
#'
is.error <- function(x) inherits(x, "try-error")

#' @title
#' Function for fitting skew-t distribution using fixed nu param
#'
#' @param vec vector of values, at least 10 data points required (much more probably better)
#' @param starting_nu Some initial value of nu >= 2
#'
#' @return log likelihood of skewed t-dist fit
#'
LRT_init <- function(vec,starting_nu){
  init_mod <- data.frame(data = vec) %>%
    selm(data ~ 1, family = "ST", fixed.param = list(nu=starting_nu), data = .)
  st_LL <- data.frame(data = vec) %>%
    selm(data ~ 1, family = "ST", data = ., start = c(coef(init_mod, param.type = 'DP'), list(nu = starting_nu))) %>%
    logLik %>%
    as.numeric()
  return(st_LL)
}

#' @title
#' LRT of skew-t compared to normal
#'
#' @description
#' Measure normality likelihood ratio statistic, as described in MacDonald et al Cell 2017
#'
#' @param vec vector of values, at least 10 data points required (much more probably better)
#'
#' @return twice log likelihood ratio of skewed t-dist fit compared to guassian fit
#' @export
#'
normLRT_test <- function(vec) {
  library(MASS); library(sn)
  min_length <- 10
  stopifnot(is.vector(vec))
  vec <- vec[!is.na(vec)]
  if (length(vec) < min_length) {return(NA)}
  st_LL <- tryCatch({
    data.frame(data = vec) %>%
      selm(data ~ 1, family = "ST", data = .) %>%
      logLik %>%
      as.numeric()
  }, error = function(e) {

    print('Default fit failed, trying set nu values')
    nu_range <- c(2, 5, 10, 25, 50, 100, 250, 500, 1000)
    i = 1; unsolved = T
    while (i <= length(nu_range) & unsolved){
      st_LL <- try(LRT_init(vec,nu_range[i]), silent=T)
      unsolved <- is.error(st_LL)
      i = i + 1
    }
    if (unsolved){
      return(NA)
    } else {
      return(st_LL)
    }
  })
  if (!is.na(st_LL)){
    n_LL <- fitdistr(vec, 'normal')$loglik
    return(2*(st_LL - n_LL))
  } else {
    return(NA)
  }
}

