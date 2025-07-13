model <- '
factor1 =~ .7*ind1 + .7*ind2 + .7*ind3 + .7*ind4
factor2 =~ .7*ind5 + .7*ind6 + .7*ind7 + .7*ind8
factor1 ~~ .3*factor2
factor1 ~~ 1*factor1
factor2 ~~ 1*factor2
'

mlr_power_bin_search <- function(model = NULL,
                                 desired_power = 0.80,
                                 alpha = 0.05,
                                 datasets = 750,
                                 left = 0,
                                 right = 1000,
                                 verbose = TRUE) {
  cl <- match.call()
  args_list <- as.list(cl)[-1]

  first_left = left
  first_right = right
  sample_range = first_right - first_left

  if(verbose == TRUE){
    number_line <- rep("-",times = 100)
    cat(paste(left,paste(number_line,collapse = ""),right,"\n",sep=" "))
  }

  power_args <- args_list[grep("^b|rx|alpha|datasets",names(args_list))]
  # Exploring the search space

  while (left <= right) {
    middle = round((left + right) / 2)

    middle_power = do.call(what = mlr_power_calc, # this is where the model power calculation function goes
                           args = as.list(c(N = middle, power_args)))
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
