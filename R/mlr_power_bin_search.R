#' Power analysis via binary search
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
#' @param desired_power The desired power for the test of the first beta coefficient.
#' @param alpha Type 1 error rate for the significance test of the first beta coefficient.
#' @param datasets The number of generated data sets used to calculate power empirically.
#' @param left This is the lower bound of the sample size search space.
#' @param right This is the upper bound of the sample size search space.
#' @param verbose Set equal to TRUE if you want to see a visualization of the search algorithm.
#' @returns A list with one element named 'N' which is the sample size needed.
#' @export
#'
#' @examples mlr_power_bin_search(b1 = 0.3, b2 = 0.1, rx1.x2 = 0.2, desired_power = 0.80, alpha = 0.05)
mlr_power_bin_search <- function(
                                b1 = 0.0,
                                b2 = 0.0,
                                b3 = 0.0,
                                b4 = 0.0,
                                b5 = 0.0,
                                rx1.x2 = 0.0,
                                rx1.x3 = 0.0,
                                rx1.x4 = 0.0,
                                rx1.x5 = 0.0,
                                rx2.x3 = 0.0,
                                rx2.x4 = 0.0,
                                rx2.x5 = 0.0,
                                rx3.x4 = 0.0,
                                rx3.x5 = 0.0,
                                rx4.x5 = 0.0,
                                desired_power = 0.80,
                                alpha = 0.05,
                                datasets = 750,
                                left = 0,
                                right = 1000,
                                verbose = TRUE){
  cl <- match.call()
  args_list <- as.list(cl)[-1]

  betas <- c(b1,b2,b3,b4,b5)
  correlations_between_preds <- c(rx1.x2,rx1.x3,rx1.x4,rx1.x5,rx2.x3,rx2.x4,rx2.x5,rx3.x4,rx3.x5,rx4.x5)

  if(b1 == 0){
    stop(call. = FALSE, "You did not specify a value for b1!")
  }
  if(any(alpha <= 0, alpha >= 1)){
    stop(call. = TRUE, "Type 1 error rate must be between 0 and 1!")
  }
  if(any(desired_power <= 0, desired_power >= 1)){
    stop(call. = TRUE, "Desired power must be between 0 and 1!")
  }
  if(any(abs(betas) >= 1)){
    stop(call. = TRUE, "Coefficient estimates must be less than 1 in magnitude!")
  }
  if(any(abs(correlations_between_preds) >= 1)){
    stop(call. = TRUE, "The correlations between predictors must be less than 1 in magnitude!")
  }

  first_left = left
  first_right = right
  sam_range = first_right - first_left

  if(verbose == TRUE){
    number_line <- rep("-",times = 100)
    cat(paste(left,paste(number_line,collapse = ""),right,"\n",sep=" "))
    }

  power_args <- args_list[grep("^b|rx|alpha|datasets",names(args_list))]
  # Exploring the search space

  while (left <= right) {
    middle = round((left + right) / 2)

    middle_power = do.call(what=mlr_power_calc, # this is where the model power calculation function goes
                           args=as.list(c(N = middle, power_args)))
    if(verbose == TRUE){
      cat(paste(rep(" ",times = length(number_line) * ((middle - first_left)/sam_range)), collapse = ""),
        "N = ", middle,
        paste(rep(" ", times = length(number_line) * (1 - ((middle - first_left)/sam_range))),collapse = ""),"\r",sep = "")
      }

    if (middle_power < desired_power) {

      left = middle + 1

    } else if (middle_power > desired_power + 0.01) {

      right = middle - 1

    } else if (middle_power > desired_power &
               middle_power < desired_power + 0.01) {

      if(verbose == TRUE){
        cat(paste(rep(" ",times = length(number_line) * ((middle - first_left)/sam_range)), collapse = ""),
          "N = ",middle,
          paste(rep(" ", times = length(number_line) * (1 - ((middle - first_left)/sam_range))),collapse = ""),"\r",sep = "")
        cat("\n\nYou need ",middle," participants for ",desired_power*100,"% power to detect a standardized beta coefficient of ",b1,".\n\n",sep = "")
      }
      return(middle)
    }
  }

  # Error if bounds of sample size search space is not large enough

  if(middle == first_left | middle == first_right){
    stop(call. = TRUE, "Adjust the bounds of your search space!")
  }

  # need this condition because otherwise the algorithm might not be able to find point of 80% power.

  if (middle == right|middle==left){
    if(verbose == TRUE){
      cat(paste(rep(" ",times = length(number_line) * ((middle - first_left)/sam_range)), collapse = ""),
        "N = ",middle,
        paste(rep(" ", times = length(number_line) * (1 - ((middle - first_left)/sam_range))),collapse = ""),"\r",sep = "")
      cat("\n\nYou need ",middle," participants for ",desired_power*100,"% power to detect a standardized beta coefficient of ",b1,".\n\n",sep = "")
    }
    return(middle)
  }
}
