#' Create correlation matrix based on standardized parameter estimates in MLR model
#'
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
#'
#' @returns A list containing a vector of means and a matrix of correlations between the variables in the model.
#' @export
#'
#' @examples mlr_params_data_gen(b1 = 0.3, b2 = 0.1, rx1.x2 = 0.2)
mlr_params_data_gen <- function(b1 = 0,
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
                                rx4.x5 = 0) {

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

  # vcov matrix creation
  mus = rep(0, 6)

  S = matrix(c(rx1.x1,rx1.x2,rx3.x1,rx1.x4,rx1.x5,rx1.y ,
               rx2.x1,rx2.x2,rx2.x3,rx2.x4,rx2.x5,rx2.y ,
               rx3.x1,rx3.x2,rx3.x3,rx3.x4,rx3.x5,rx3.y ,
               rx4.x1,rx4.x2,rx4.x3,rx4.x4,rx4.x5,rx4.y ,
               rx5.x1,rx5.x2,rx5.x3,rx5.x4,rx5.x5,rx5.y ,
               rx1.y ,rx2.y ,rx3.y ,rx4.y ,rx5.y ,ry.y ),
             nrow = 6,
             byrow = TRUE)

  return(list(mus = mus, S = S))
}
