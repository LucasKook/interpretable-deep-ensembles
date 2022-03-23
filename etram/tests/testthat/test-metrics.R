
context("metrics")

set.seed(26448)

cdf1 <- data.frame(matrix(c(0.1, 0.2, 1,
                            0.5, 0.7, 1,
                            0.4, 0.7, 1),
                          nrow = 3, byrow = TRUE))
cdf2 <- data.frame(matrix(c(0.5, 0.7, 1,
                            0.2, 0.8, 1,
                            0.8, 0.9, 1),
                          nrow = 3, byrow = TRUE))
y_true <- data.frame(matrix(c(0, 0, 1,
                              1, 0, 0,
                              1, 0, 0),
                            nrow = 3, byrow = TRUE))
lys_cdf = list(cdf1, cdf2)

test_that("accuracy is correct", {
  expect_equal(get_acc(cdf = cdf1, y_true = y_true), 1)
  expect_equal(get_acc(cdf = cdf2, y_true = y_true), 1/3)
  expect_equal(get_avg_acc(lys_cdf = lys_cdf, y_true = y_true),
               mean(c(get_acc(cdf = cdf1, y_true = y_true),
                      get_acc(cdf = cdf2, y_true = y_true))))
})

test_that("rps is correct", {
  expect_equal(get_rps(cdf = cdf1, y_true = y_true),
               ((0.1^2 + 0.2^2 + 0.5^2 + 0.3^2 + 0.6^2 + 0.3^2) / 2) / 3)
  expect_equal(get_rps(cdf = cdf2, y_true = y_true),
               ((0.5^2 + 0.7^2 + 0.8^2 + 0.2^2 + 0.2^2 + 0.1^2) / 2) / 3)
  expect_equal(get_avg_rps(lys_cdf = lys_cdf, y_true = y_true),
               mean(c(get_rps(cdf = cdf1, y_true = y_true),
                      get_rps(cdf = cdf2, y_true = y_true))))
})

test_that("nll is correct", {
  expect_equal(get_nll(cdf = cdf1, y_true = y_true),
               -mean(log(c(0.8, 0.5, 0.4))))
  expect_equal(get_nll(cdf = cdf2, y_true = y_true),
               -mean(log(c(0.3, 0.2, 0.8))))
  expect_equal(get_avg_nll(lys_cdf = lys_cdf, y_true = y_true),
               mean(c(get_nll(cdf = cdf1, y_true = y_true),
                      get_nll(cdf = cdf2, y_true = y_true))))
})

test_that("binnll is correct", {
  expect_equal(get_binnll(cdf = cdf1, y_true = y_true, cutoff = 2),
               -mean(log(c(0.8, 0.7, 0.7))))
  expect_equal(get_binnll(cdf = cdf2, y_true = y_true, cutoff = 2),
               -mean(log(c(0.3, 0.8, 0.9))))
  expect_equal(get_avg_binnll(lys_cdf = lys_cdf, y_true = y_true, cutoff = 2),
               mean(c(get_binnll(cdf = cdf1, y_true = y_true, cutoff = 2),
                      get_binnll(cdf = cdf2, y_true = y_true, cutoff = 2))))
})

p1 <- runif(100)
p2 <- runif(100)
cdf1 <- data.frame(matrix(c(p1, rep(1, 100)),
                         nrow = 100, byrow = FALSE))
cdf2 <- data.frame(matrix(c(p2, rep(1, 100)),
                          nrow = 100, byrow = FALSE))
lys_cdf <- list(cdf1, cdf2)
y_true <- data.frame(t(replicate(100, sample(c(0, 1)))))

test_that("auc is correct", {
  expect_equal(get_auc(cdf = cdf1, y_true = y_true, cutoff = 1),
               c(pROC::auc(y_true[, 1], p1, levels = c(0, 1), direction = "<")),
               tolerance = 1e-4)
  expect_equal(get_auc(cdf = cdf2, y_true = y_true, cutoff = 1),
               c(pROC::auc(y_true[, 1], p2, levels = c(0, 1), direction = "<")),
               tolerance = 1e-4)
  expect_equal(get_avg_auc(lys_cdf = lys_cdf, y_true = y_true, cutoff = 1),
               mean(c(get_auc(cdf = cdf1, y_true = y_true, cutoff = 1),
                      get_auc(cdf = cdf2, y_true = y_true, cutoff = 1))))
})

test_that("cal is correct", {
  expect_equal(get_cal_perclass(cdf1, y_true)$cint,
               c(as.numeric(coef(glm(y_true[, 1] ~ offset(qlogis(cdf1[, 1])),
                                     family = "binomial"))[1]),
                 as.numeric(coef(glm(y_true[, 2] ~ offset(qlogis(1 - cdf1[, 1])),
                                     family = "binomial"))[1])))
  expect_equal(get_cal_perclass(cdf1, y_true)$cslope,
               c(as.numeric(coef(glm(y_true[, 1] ~ qlogis(cdf1[, 1]),
                                     family = "binomial"))[2]),
                 as.numeric(coef(glm(y_true[, 2] ~ qlogis(1 - cdf1[, 1]),
                                     family = "binomial"))[2])))
  expect_equal(get_cal(cdf1, y_true)$cint,
               mean(get_cal_perclass(cdf1, y_true)$cint))
  expect_equal(get_avg_cal(lys_cdf, y_true)$cint,
               mean(c(get_cal(cdf1, y_true)$cint, get_cal(cdf2, y_true)$cint)))
})
