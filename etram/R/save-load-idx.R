
#' Save row indices of random splits
#' @description Saves row indices used for training, validating and testing for each split as a data frame.
#' Ensures that all outcome classes are present in the training set.
#' @param nsplits number of splits.
#' @param prptest proportion of observations used for testing.
#' @param prpval proportion of observations used for validation.
#' @param tab_dat data frame containing at least the response variable.
#' @param fml model formula (only used to identify the response variable).
#' @export
save_ridx <- function(nsplts, prptest, prpval, tab_dat, fml, out_dir, fname) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir)
  }
  n <- nrow(tab_dat)
  y <- model.frame(fml, data = tab_dat)[, 1L]
  K <- length(unique(y))
  ridx <- seq_len(n)
  ret <- data.frame(idx = numeric(),
                    type = character(),
                    spl = numeric())
  for (spl in seq_len(nsplts)) {
    not_all_classes <- TRUE
    while (not_all_classes) {
      ridx_test <- sample(ridx, ceiling(prptest * n), replace = FALSE)
      ridx_val <- sample(ridx[!(ridx %in% ridx_test)], ceiling(prpval * n), replace = FALSE)
      ridx_train <- ridx[!(ridx %in% c(ridx_test, ridx_val))]
      tmp <- data.frame(idx = c(ridx_test, ridx_val, ridx_train),
                        type = c(rep("test", length(ridx_test)),
                                 rep("val", length(ridx_val)),
                                 rep("train", length(ridx_train))),
                        spl = spl)
      if (length(unique(y[ridx_train])) == K) {
        not_all_classes <- FALSE
      }
    }
    ret <- rbind(ret, tmp)
  }
  write.csv(ret, file = paste0(out_dir, fname, "_ridx.csv"))
}

#' Load row indices of random splits
#' @export
get_ridx <- function(in_dir, fname) {
  read.csv(file = paste0(in_dir, fname, "_ridx.csv"))[, -1L]
}
