### Movies experiment
### BS/LK 2023

# Dependencies ------------------------------------------------------------

library("deeptrafo")
library("tidyverse")
library("ggpubr")

# Params ------------------------------------------------------------------

bpath <- "scratch/movie-exp-JMLR"
nr_words <- 10000 # Number of words used to describe the text
embedding_size <- 99 # Size of text embedding
repl <- 1 # 3554 movies

do_fit <- FALSE

# Prepare data ------------------------------------------------------------

dl <- readRDS(file.path(bpath, "data_splitted.RDS"))

mtrain <- tibble(
  log_revenue = log10(1 + dl[[repl]]$train$revenue),
  log_popularity = log10(dl[[repl]]$train$popularity),
  texts = dl[[repl]]$train$texts
) %>% filter(log_revenue > 4)

test <- tibble(
  log_revenue = log10(1 + dl[[repl]]$test$revenue),
  log_popularity = log10(dl[[repl]]$test$popularity),
  texts = dl[[repl]]$test$texts
) %>% filter(log_revenue > 4)

### split train in train and valid for ens-weight tuning
set.seed(123)
train_size <- round(nrow(mtrain) * 0.95)
train_indices <- sample.int(nrow(mtrain), train_size, replace = FALSE)
valid_indices <- setdiff(seq_len(nrow(mtrain)), train_indices)
train <- mtrain[train_indices, ]
valid <- mtrain[valid_indices, ]

# prepare grid for plots
no_grid <- 150
no_test <- nrow(test)
grid <- seq(min(train$log_revenue), max(train$log_revenue), length.out = no_grid)

### prepare test w response grid
tnr <- test %>%
  rename(response = log_revenue) %>%
  mutate(idx = 1:nrow(.)) %>%
  nest_by(idx) %>%
  mutate(log_revenue = list(grid)) %>%
  unnest(c(data, log_revenue))

# Define neural network architecture
nn <- \(x) x %>%
  layer_embedding(input_dim = nr_words, output_dim = embedding_size) %>%
  layer_lstm(units = 50L, return_sequences = TRUE) %>%
  layer_lstm(units = 50L, return_sequences = FALSE) %>%
  layer_dropout(rate = 0.1) %>%
  layer_dense(25L) %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(5L, name = 'penultimate') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(1L)

## fit ensemble semi Model ##########################

# ens <- ensemble(m_semi, n_ensemble = 3, epochs = 0) # test

no_epochs <- 120 # 15 # 50
no_ens <- 5

my_seed <- 123867
(fileName <- file.path(bpath, paste0("weights_ensemble_seed_", my_seed,".RDS")))
seed_offset <- my_seed

if (do_fit) {
  ens <- trafoensemble(
    log_revenue ~ 0 + log_popularity + nn(texts), data = train,
    latent_distr = "logistic",
    list_of_deep_models = list(nn = nn),
    validation_split = 0.1,
    callbacks = list(callback_early_stopping(patience = 7, restore_best_weights = TRUE)),
    epochs = no_epochs, n_ensemble = no_ens,
    tf_seeds = seq(seed_offset, (seed_offset + no_ens - 1), by = 1),
    seed = seq(seed_offset, (seed_offset + no_ens - 1), by = 1), verbose = TRUE
  )

  saveRDS(ens$ensemble_results, fileName) # save as rds
} else {
  # step 1: built model w/o training: ens <- ... # with epochs = 0
  ens <- trafoensemble(
    log_revenue ~ 0 + log_popularity + nn(texts), data = train,
    latent_distr = "logistic", list_of_deep_models = list(nn = nn),
    validation_split = 0.1, callbacks = list(
      callback_early_stopping(patience = 7, restore_best_weights = TRUE)),
    epochs = 0, n_ensemble = no_ens, tf_seeds = seq(
      seed_offset, (seed_offset + no_ens - 1), by = 1), seed = seq(
        seed_offset, (seed_offset + no_ens - 1), by = 1), verbose = TRUE
  )

  # step 2: load weights from saved RDS
  ens$ensemble_results <- readRDS(fileName)
}

logLik(ens, newdata = test, convert_fun = mean, batch_size = 64)
tuned <- weighted_logLik(ens, newdata = valid) # val dat for tuning weights
weighted_logLik(ens, weights = tuned$weights, newdata = test)

# Vis ---------------------------------------------------------------------

fitted <- simplify2array(fitted(ens, type = "trafo", newdata = tnr, batch_size = 1e4))
fitens <- apply(fitted, 1:2, mean)

tnr$trafo <- fitens[, 2] + fitens[, 1]
tnr$dens <- dlogis(tnr$trafo) * fitens[, 4]
tnr$cdf <- plogis(tnr$trafo)

mlwd <- 1
kyhigh <- unique(tnr$idx[tnr$response > 8])[6]
kylow <- unique(tnr$idx[tnr$response < 6])[6]

col_low = "darkblue"
col_high = "red"

col_m_low = "skyblue"
col_m_high = "orange"

dhigh <- filter(tnr, idx == kyhigh)
dlow <- filter(tnr, idx == kylow)

dlow$mempdf <- do.call("cbind", predict(ens, newdata = dlow, type = "pdf"))
dhigh$mempdf <- do.call("cbind", predict(ens, newdata = dhigh, type = "pdf"))

dlow$memcdf <- do.call("cbind", predict(ens, type = "cdf", newdata = dlow))
dhigh$memcdf <- do.call("cbind", predict(ens, type = "cdf", newdata = dhigh))

p1 <- ggplot(test, aes(x = log_popularity, y = log_revenue)) +
  geom_point(alpha = 0.3) +
  geom_point(aes(y = response), data = dhigh, color = col_high) +
  geom_point(aes(y = response), data = dlow, color = col_low) +
  geom_hline(yintercept = dhigh$response[1], color = col_high, linetype = 2) +
  geom_hline(yintercept = dlow$response[1], color = col_low, linetype = 2) +
  theme_bw() + labs(x = expression(Popularity~score~(log[10]*'-'*scale)),
                    y = expression(Movie~revenue~(log[10]*'-'*scale)))

p2 <- ggplot(dhigh, aes(x = log_revenue, y = dens)) +
  geom_line(aes(y = mempdf[, 1], linetype = "member"), col = col_m_high) +
  geom_line(aes(y = mempdf[, 2], linetype = "member"), col = col_m_high) +
  geom_line(aes(y = mempdf[, 3], linetype = "member"), col = col_m_high) +
  geom_line(aes(y = mempdf[, 4], linetype = "member"), col = col_m_high) +
  geom_line(aes(y = mempdf[, 5], linetype = "member"), col = col_m_high) +
  geom_line(lwd = mlwd, aes(color = "high", linetype = "ensemble")) +
  geom_line(aes(y = mempdf[, 1], linetype = "member"), data = dlow, col = col_m_low) +
  geom_line(aes(y = mempdf[, 2], linetype = "member"), data = dlow, col = col_m_low) +
  geom_line(aes(y = mempdf[, 3], linetype = "member"), data = dlow, col = col_m_low) +
  geom_line(aes(y = mempdf[, 4], linetype = "member"), data = dlow, col = col_m_low) +
  geom_line(aes(y = mempdf[, 5], linetype = "member"), data = dlow, col = col_m_low) +
  geom_line(lwd = mlwd, data = dlow, aes(color = "low", linetype = "ensemble")) +
  scale_color_manual(values = c(col_high, col_low)) +
  theme_bw() +
  labs(y = expression(Conditional~density~of~Y~"|"~X==x),
       x = expression(Movie~revenue~(log[10]*'-'*scale)),
       color = "Observation", linetype = element_blank()) +
  geom_vline(xintercept = c(dlow$response[1], dhigh$response[1]),
             col = c(col_low, col_high), linetype = 2)

dhigh$avg <- apply(-log(dhigh$mempdf), 1, mean)
dlow$avg <- apply(-log(dlow$mempdf), 1, mean)

p4 <- ggplot(dhigh, aes(x = log_revenue, y = -log(dens))) +
  geom_line(aes(y = -log(mempdf[, 1])), linetype = 3, col = col_m_high) +
  geom_line(aes(y = -log(mempdf[, 2])), linetype = 3, col = col_m_high) +
  geom_line(aes(y = -log(mempdf[, 3])), linetype = 3, col = col_m_high) +
  geom_line(aes(y = -log(mempdf[, 4])), linetype = 3, col = col_m_high) +
  geom_line(aes(y = -log(mempdf[, 5])), linetype = 3, col = col_m_high) +
  # geom_line(aes(y = avg), lwd = 1, col = "gray60", data = dhigh) +
  geom_line(lwd = mlwd, col = col_high) +
  geom_line(aes(y = -log(mempdf[, 1])), linetype = 3, data = dlow, col = col_m_low) +
  geom_line(aes(y = -log(mempdf[, 2])), linetype = 3, data = dlow, col = col_m_low) +
  geom_line(aes(y = -log(mempdf[, 3])), linetype = 3, data = dlow, col = col_m_low) +
  geom_line(aes(y = -log(mempdf[, 4])), linetype = 3, data = dlow, col = col_m_low) +
  geom_line(aes(y = -log(mempdf[, 5])), linetype = 3, data = dlow, col = col_m_low) +
  # geom_line(aes(y = avg), lwd = 1, data = dlow, col = "gray60") +
  geom_line(lwd = mlwd, data = dlow, col = col_low) +
  theme_bw() +
  labs(y = expression(Negative~log*'-'*likelihood),
       x = expression(Movie~revenue~(log[10]*'-'*scale))) +
  geom_vline(xintercept = c(dlow$response[1], dhigh$response[1]),
             col = c(col_low, col_high), linetype = 2)

cl <- apply(
  tr <- do.call("cbind", predict(ens, newdata = test, type = "cdf")), 1, mean
)

p3 <- data.frame(cll = qlogis(cl), trr = qlogis(tr)) %>%
  pivot_longer(trr.1:trr.5, names_to = "member", values_to = "trafo") %>%
  ggplot(aes(x = cll, y = trafo, color = member)) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(y = expression(Transformation~h[m](y~"|"~x)),
       x = expression(Classical~ensemble~{F^{-1}}[Z]({bar(F)^c}[M](y~"|"~x))))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE)
ggsave(file.path(bpath, "figure4.pdf"), height = 6, width = 7)
