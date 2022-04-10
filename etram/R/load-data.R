
#' Load data
#' @export
load_data <- function(which = c("mnist", "stroke", "utkface", "melanoma"), path = NULL, im_path = NULL) {
  which <- match.arg(which)

  if (which == "mnist") {
    mnist <- dataset_mnist()
    c(c(x_train, y_train), c(x_test, y_test)) %<-% mnist
    im <- array_reshape(x_train, c(60000, 28, 28, 1))
    im <- im / 255
    tab_dat <- data.frame(y = y_train)
    tab_dat$y <- factor(tab_dat$y, ordered = TRUE)
    ret <- list(tab_dat = tab_dat,
                im = im)

  } else if (which == "stroke") {
    im <- h5read(file = im_path, "/X")
    im <- array(aperm(im, 4:1), dim = c(dim(im)[4:1], 1))
    y <- h5read(file = im_path, "/Y_pat")
    pat <- data.frame(p_id = h5read(file = im_path, "/pat"))
    tab_dat <- read_csv(path, na = c("NA")) %>%
      left_join(pat, .) %>%
      mutate(mrs3 = ordered(mrs3, levels = 0:6),
             mrs_before = factor(mrs_before, levels = unique(na.omit(mrs_before)),
                                 labels = 0:4),
             mrs3_unfavorable = factor(ifelse(mrs3 <= "2", 0, 1), ordered = TRUE))
    tab_dat <- tab_dat[complete.cases(tab_dat), ]
    im <- ontram:::.batch_subset(im, as.numeric(row.names(tab_dat)), dim(im))
    ret <- list(tab_dat = tab_dat,
                im = im)

  } else if (which == "utkface") {
    train <- .read_utkface("train", path = path)
    valid <- .read_utkface("valid", path = path)
    test <- .read_utkface("test", path = path)
    im_train <- train[[1]]
    im_val <- valid[[1]]
    im_test <- test[[1]]
    tab_train <- train[[2]]
    tab_val <- valid[[2]]
    tab_test <- test[[2]]
    colnames(tab_train) <- colnames(tab_val) <- colnames(tab_test) <- c("age_group", "gender", "race",
                                                                        "x_1", "x_2", "x_3", "x_4", "x_5",
                                                                        "x_6", "x_7", "x_8", "x_9", "x_10")
    tab_dat <- rbind(tab_train, tab_val, tab_test)
    tab_dat$age_group <- factor(tab_dat$age_group, levels = 0:6, ordered = TRUE)
    im <- abind::abind(im_train, im_val, im_test, along = 1)
    ret <- list(tab_dat = tab_dat,
                im = im)

  } else if (which == "melanoma") {
    tab_dat <- read_csv(path) %>%
      mutate(target = factor(target))
    tab_dat <- tab_dat[!is.na(tab_dat$age_approx), ] # exclude missings
    tab_dat$target <- ordered(tab_dat$target, levels = 0:1)
    # standardization
    m <- mean(tab_dat$age_approx)
    s <- sd(tab_dat$age_approx)
    tab_dat$age_s = (tab_dat$age_approx - m)/(s + 1e-12)
    tab_dat <- as.data.frame(tab_dat)
    image_size <- 128L
    im <- array(dim = c(nrow(tab_dat), image_size, image_size, 3L))
    for (i in seq_len(nrow(tab_dat))) {
      print(i)
      name <- paste0(im_path, tab_dat$image_name[i], '.jpg')
      if (!file.exists(name)){
        print(paste0('File ', name, ' does not exist'))
      }
      img <- image_load(name, target_size = c(image_size, image_size))
      img <- image_to_array(img)
      im[i,,,] <- img/255
    }
    ret <- list(tab_dat = tab_dat,
                im = im)
  }
  return(ret)
}

.read_utkface <- function(mode, path) {
  dset <- paste0("_", mode)
  inm <- paste0("X", dset) # images
  nms <- paste0(c("age_group", "gender", "race",
                  paste0("x_", 1:10)), dset) # clinical data and response
  im <- h5read(path, inm)
  im <- aperm(im, 4:1)
  dat <- sapply(nms, function(nm) {h5read(path, nm)})
  dat <- data.frame(dat)
  ret <- list(im, dat)
  return(ret)
}
