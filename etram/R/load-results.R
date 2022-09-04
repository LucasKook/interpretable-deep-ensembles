
#' Load CDFs of all splits and ensemble members
#' @export
list_cdfs <- function(in_dir, fname, splits, ensembles = NULL, type = c("test", "val", "train")) {
  type <- match.arg(type)
  lys <- list()
  for (spl in seq_len(splits)) {
    eval(parse(text = paste0("lys", spl, "<- list()")))
    if (!is.null(ensembles)) {
      for (ens in seq_len(ensembles)) {
        assign(paste0("cdf", ens),
               read.csv(file.path(in_dir, paste0(fname, "_cdf", type, "_spl",
                                                 spl, "_ens", ens, ".csv")))[, -c(1L, 2L), drop = FALSE])
        eval(parse(text = paste0("lys", spl, "[[", ens, "]] <- cdf", ens)))
      }
    } else { # si, sils
      assign(paste0("cdf", spl),
             read.csv(file.path(in_dir, paste0(fname, "_cdf", type, "_spl", spl, ".csv")))[, -c(1L, 2L), drop = FALSE])
      eval(parse(text = paste0("lys", spl, "[[1]] <- cdf", spl)))
    }
    eval(parse(text = paste0("lys", "[[", spl, "]] <- lys", spl)))
  }
  return(lys)
}

#' Load predictions on the scale of \code{h} of all splits and ensemble members
#' @export
list_trafos <- function(in_dir, fname, splits, ensembles = NULL, type = c("test", "val", "train")) {
  type <- match.arg(type)
  lys <- list()
  for (spl in seq_len(splits)) {
    eval(parse(text = paste0("lys", spl, "<- list()")))
    if (!is.null(ensembles)) {
      for (ens in seq_len(ensembles)) {
        assign(paste0("trafo", ens),
               read.csv(file.path(in_dir, paste0(fname, "_trafo", type, "_spl",
                                                 spl, "_ens", ens, ".csv")))[, -c(1L), drop = FALSE])
        eval(parse(text = paste0("lys", spl, "[[", ens, "]] <- trafo", ens)))
      }
    } else { # si, sils
      assign(paste0("trafo", spl),
             read.csv(file.path(in_dir, file.path(fname, "_trafo", type, "_spl",
                                                  spl, ".csv")))[, -c(1L), drop = FALSE])
      eval(parse(text = paste0("lys", spl, "[[1]] <- trafo", spl)))
    }
    eval(parse(text = paste0("lys", "[[", spl, "]] <- lys", spl)))
  }
  return(lys)
}

#' Load raw terms on the scale of \code{h} of all splits and ensemble members
#' @export
list_terms <- function(in_dir, fname, splits, ensembles = NULL, type = c("test", "val", "train")) {
  type <- match.arg(type)
  lys <- list()
  for (spl in seq_len(splits)) {
    eval(parse(text = paste0("lys", spl, "<- list()")))
    if (!is.null(ensembles)) {
      for (ens in seq_len(ensembles)) {
        assign(paste0("raw", ens),
               read.csv(file.path(in_dir, paste0(fname, "_raw", type, "_spl",
                                                 spl, "_ens", ens, ".csv")))[, -c(1L), drop = FALSE])
        eval(parse(text = paste0("lys", spl, "[[", ens, "]] <- raw", ens)))
      }
    } else { # si, sils
      assign(paste0("raw", spl),
             read.csv(file.path(in_dir, paste0(fname, "_raw", type, "_spl",
                                               spl, ".csv")))[, -c(1L), drop = FALSE])
      eval(parse(text = paste0("lys", spl, "[[1]] <- raw", spl)))
    }
    eval(parse(text = paste0("lys", "[[", spl, "]] <- lys", spl)))
  }
  return(lys)
}

#' Load log-odds ratios of all splits and ensemble members
#' @export
list_lors <- function(in_dir, fname, splits, ensembles = NULL) {
  lys <- list()
  for (spl in seq_len(splits)) {
    eval(parse(text = paste0("lys", spl, "<- list()")))
    if (!is.null(ensembles)) {
      for (ens in seq_len(ensembles)) {
        assign(paste0("lor", ens),
               read.csv(file.path(in_dir, paste0(fname, "_lor_spl", spl, "_ens", ens, ".csv"))))
        eval(parse(text = paste0("lys", spl, "[[", ens, "]] <- lor", ens)))
      }
    } else { # si, sils
      assign(paste0("lor", spl),
             read.csv(file.path(in_dir, paste0(fname, "_lor_spl", spl, ".csv"))))
      eval(parse(text = paste0("lys", spl, "[[1]] <- lor", spl)))
    }
    eval(parse(text = paste0("lys", "[[", spl, "]] <- lys", spl)))
  }
  return(lys)
}

#' Load histories of all splits and ensemble members
#' @export
load_hists <- function(in_dir, fname, splits, ensembles) {
  lys <- list()
  for (spl in seq_len(splits)) {
    eval(parse(text = paste0("lys", spl, "<- list()")))
    for (ens in seq_len(ensembles)) {
      assign(paste0("hist", ens),
             read.csv(file = file.path(in_dir, paste0(fname, "_hist_spl", spl,
                                                      "_ens", ens, ".csv"))))
      eval(parse(text = paste0("hist", ens, "$Split <- ", spl)))
      eval(parse(text = paste0("hist", ens, "$Ens <- ", ens)))
      eval(parse(text = paste0("lys", spl, "[[", ens, "]]", "<- hist", ens)))
    }
    eval(parse(text = paste0("df", spl, "<- do.call(rbind, lys", spl, ")")))
    eval(parse(text = paste0("lys", "[[", spl, "]] <- df", spl)))
  }
  df <- do.call(rbind, lys)
  df$Type <- as.factor(df$Type)
  df$Split <- as.factor(df$Split)
  df$Ens <- as.factor(df$Ens)
  return(df)
}
