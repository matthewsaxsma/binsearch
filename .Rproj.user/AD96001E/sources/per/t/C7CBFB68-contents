#' Power of test for first beta coefficient
#'
#' @param N The sample size.
#' @param mus Vector of means.
#' @param S Correlation matrix of variables.
#' @param alpha Type 1 error rate for the significance test of the first beta ceofficient.
#' @param datasets The number of generated datasets used to calculate power empirically.
#'
#' @returns A numeric vector of length 1 which contains the power of the test.
#' @export
#' @import MASS
#' @examples
#' mus = c(0,0,0)
#' S = matrix(c(1.0,0.2,0.1,
#'              0.2,1.0,0.3,
#'              0.1,0.3,1.0), nrow = 3, byrow = TRUE)
#' MLR_power_calc(N = 150, mus = mus, S = S, alpha = 0.05, datasets = 1000)
MLR_power_calc <- function(N = NULL,
                           mus = NULL,
                           S = NULL,
                           alpha = NULL,
                           datasets = NULL) {
  pval = numeric(datasets)
  modelformula = as.formula(paste("y ~", paste(
    "x", 1:(ncol(S)-1), sep = "", collapse = " + "
  )))
  for (i in 1:datasets) {
    # Data generation
    temp = data.frame(MASS::mvrnorm(N,mu = mus,Sigma = S,empirical = FALSE))
    names(temp) = c(paste("x", seq(1, (ncol(S)-1)), sep = ""), "y")
    temp = as.data.frame(scale(temp, center = TRUE, scale = TRUE))
    # Need to solve (X'X)^-1X'Y
    fit = lm(modelformula, data = temp)
    # b1 coefficient p-value
    pval[i] = summary(fit)$coefficients["x1", "Pr(>|t|)"]
  }
  power = mean(pval < alpha)
  return(power)
}
