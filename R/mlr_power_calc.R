#' Power calculation for coefficient of multiple linear regression model
#'
#' @param N The sample size
#' @param b1 The coefficient you want to be powered to detect.
#' @param b2 A second beta coefficient.
#' @param b3 A third beta coefficient.
#' @param b4 A fourth beta coefficient.
#' @param b5 The fifth beta coefficient.
#' @param rx1.x2 The correlation between predictors 1 and 2.
#' @param rx1.x3 The correlation between predictors 1 and 3.
#' @param rx1.x4 The correlation between predictors 1 and 4.
#' @param rx1.x5 The correlation between predictors 1 and 5.
#' @param rx2.x3 The correlation between predictors 2 and 3.
#' @param rx2.x4 The correlation between predictors 2 and 4.
#' @param rx2.x5 The correlation between predictors 2 and 5.
#' @param rx3.x4 The correlation between predictors 3 and 4.
#' @param rx3.x5 The correlation between predictors 3 and 5.
#' @param rx4.x5 The correlation between predictors 4 and 5.
#' @param alpha Type 1 error rate for the significance test of the first beta ceofficient.
#' @param datasets The number of generated datasets used to calculate power empirically.
#'
#' @return Returns the power of the test
#' @export
#'
#' @examples mlr_power_calc(N = 600, b1 = 0.1, b2 = 0.3, rx1.x2 = 0.25)


mlr_power_calc <- function(N = 0,
                           b1 = 0,
                           b2 = 0,
                           b3 = 0,
                           b4 = 0,
                           b5 = 0,
                           rx1.x2 = 0,
                           rx1.x3 = 0,
                           rx1.x4 = 0,
                           rx1.x5 = 0,
                           rx2.x3 = 0,
                           rx2.x4 = 0,
                           rx2.x5 = 0,
                           rx3.x4 = 0,
                           rx3.x5 = 0,
                           rx4.x5 = 0,
                           alpha = 0.05,
                           datasets = 750){
  betas = c(b1,b2,b3,b4,b5)
  correlations_between_preds <- c(rx1.x2,rx1.x3,rx1.x4,rx1.x5,rx2.x3,rx2.x4,rx2.x5,rx3.x4,rx3.x5,rx4.x5)

  npred = sum(!betas==0) # soft count of predictors in model. however, they could still specify non-zero correlations between predictors where one has no effect on the outcome
  # set the variances to be 1
  rx1.x1 = 1
  rx2.x2 = 1
  rx3.x3 = 1
  rx4.x4 = 1
  rx5.x5 = 1
  ry.y = 1

  # switch variables for easy matrix creating
  rx2.x1 = rx1.x2
  rx3.x1 = rx1.x3
  rx4.x1 = rx1.x4
  rx5.x1 = rx1.x5
  rx3.x2 = rx2.x3
  rx4.x2 = rx2.x4
  rx5.x2 = rx2.x5
  rx4.x3 = rx3.x4
  rx5.x3 = rx3.x5
  rx5.x4 = rx4.x5

  # correlations between predictors and outcome
  rx1.y = b1 + rx1.x2 * b2 + rx1.x3 * b3 + rx1.x4 * b4 + rx1.x5 * b5
  rx2.y = b2 + rx1.x2 * b1 + rx2.x3 * b3 + rx2.x4 * b4 + rx2.x5 * b5
  rx3.y = b3 + rx1.x3 * b1 + rx2.x3 * b2 + rx1.x4 * b4 + rx1.x5 * b5
  rx4.y = b4 + rx1.x4 * b1 + rx2.x4 * b2 + rx3.x4 * b3 + rx4.x5 * b5
  rx5.y = b5 + rx1.x5 * b1 + rx2.x5 * b2 + rx3.x5 * b3 + rx4.x5 * b4

  # variance-covariance matrix creation

  S = matrix(c(rx1.x1,rx1.x2,rx3.x1,rx1.x4,rx1.x5,rx1.y ,
               rx2.x1,rx2.x2,rx2.x3,rx2.x4,rx2.x5,rx2.y ,
               rx3.x1,rx3.x2,rx3.x3,rx3.x4,rx3.x5,rx3.y ,
               rx4.x1,rx4.x2,rx4.x3,rx4.x4,rx4.x5,rx4.y ,
               rx5.x1,rx5.x2,rx5.x3,rx5.x4,rx5.x5,rx5.y ,
               rx1.y ,rx2.y ,rx3.y ,rx4.y ,rx5.y ,ry.y ),
             nrow = 6,
             byrow = TRUE)

  mus = rep(0, ncol(S))
  pval = numeric(datasets)
  modelformula = as.formula(paste("y ~", paste("x", 1:npred, sep = "", collapse = " + ")))
  temp = data.frame(matrix(NA,nrow = N,ncol = ncol(S)))
  names(temp) = c(paste("x", seq(1,ncol(S)-1), sep = ""), "y")

  for (i in 1:datasets) {
    # Data generation
    temp[,1:6] = MASS::mvrnorm(N,mu = mus,Sigma = S,empirical = FALSE)
    temp = scale(temp,center = TRUE,scale=TRUE)
    temp = as.data.frame(temp)
    # Need to solve (X'X)^-1X'Y
    fit = lm(modelformula, data = temp)
    # b1 coefficient p-value
    pval[i] = summary(fit)$coefficients["x1", "Pr(>|t|)"]
  }
  power = mean(pval < alpha)
  return(power)
}
