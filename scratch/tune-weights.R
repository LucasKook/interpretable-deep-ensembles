# Tuning ensemble weights in deeptrafo
# LK Apr 2023

set.seed(1964)

# Deps --------------------------------------------------------------------

library("deeptrafo")

# Data --------------------------------------------------------------------

dgp <- function(n = 1e2) {
  x = runif(n)
  y = sin(x * pi * 3) + rnorm(n, sd = 0.1)
  data.frame(y = y, x = x)
}

weighted_logLik <- function(
    object,
    weights = NULL,
    newdata = NULL,
    convert_fun = function(x, ...) mean(x, ...),
    batch_size = NULL,
    ...
) {

  indiv <- deeptrafo:::.call_for_all_members(
    object, deeptrafo:::logLik.deeptrafo, newdata = newdata, y = y,
    convert_fun = convert_fun, ... = ...
  )

  fitt <- fitted(object, newdata = newdata, batch_size = NULL)

  if (is.null(newdata)) {
    y <- object$init_params$y
  } else {
    y <- deeptrafo:::response(newdata[[object$init_params$response_varname]])
  }

  if (is.null(weights)) {
    obj <- function(weights) {
      y_pred <- apply(simplify2array(fitt), 1:2, weighted.mean, w = weights)
      convert_fun(object$model$loss(y, y_pred)$numpy())
    }
    opt <- optim(rep(1, length(fitt)), obj, lower = .Machine$double.eps,
                 upper = 1 - .Machine$double.eps, method = "L-BFGS-B")
    ensemble_loss <- opt$value
    weights <- opt$par
    weights <- weights / sum(weights)
  } else {
    y_pred <- apply(simplify2array(fitt), 1:2, weighted.mean, w = weights)
    ensemble_loss <- convert_fun(object$model$loss(y, y_pred)$numpy())
  }

  list(members = unlist(indiv),
       mean = mean(unlist(indiv)),
       ensemble = ensemble_loss,
       weights = weights)

}

### Generate train, validation, test data
train_data <- dgp()
validation_data <- dgp()
test_data <- dgp()

### Train ensemble with early stopping (separate validation split)
ens <- trafoensemble(
  y ~ x, data = train_data, epochs = 1e3, n_ensemble = 10, verbose = TRUE,
  validation_split = 0.1, optimizer = optimizer_adam(learning_rate = 0.1),
  callbacks = list(callback_early_stopping(patience = 50))
)

### Compute the optimal weights on the validation data (weights = NULL; default)
tuned <- weighted_logLik(ens, newdata = validation_data)

### Use optimal weights and make test predictions (weights = tuned$weights)
weighted_logLik(ens, weights = tuned$weights, newdata = test_data)
