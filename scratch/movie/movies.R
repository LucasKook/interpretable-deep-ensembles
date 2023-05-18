### Movies experiment
### BS/LK 2023

# Dependencies ------------------------------------------------------------

library("deeptrafo")
library("tidyverse")
library("ggpubr")

# Params ------------------------------------------------------------------

bpath <- "scratch/movie"
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

mlwd <- 0.8
kyhigh <- unique(tnr$idx[tnr$response > 8])[6]
kylow <- unique(tnr$idx[tnr$response < 6])[6]

col2 <- viridis::viridis(2, begin = 0.1, end = 0.9)
col_low <- col2[1]
col_high <- col2[2]

col_m_low <- "skyblue"
col_m_high <- "orange"

dhigh <- filter(tnr, idx == kyhigh)
dlow <- filter(tnr, idx == kylow)

dlow$mempdf <- do.call("cbind", predict(ens, newdata = dlow, type = "pdf"))
dhigh$mempdf <- do.call("cbind", predict(ens, newdata = dhigh, type = "pdf"))
dlow$enslin <- rowMeans(dlow$mempdf)
dhigh$enslin <- rowMeans(dhigh$mempdf)

dlow$memcdf <- do.call("cbind", predict(ens, type = "cdf", newdata = dlow))
dhigh$memcdf <- do.call("cbind", predict(ens, type = "cdf", newdata = dhigh))
dlow$enslincdf <- rowMeans(dlow$memcdf)
dhigh$enslincdf <- rowMeans(dhigh$memcdf)

dhigh$favg <- exp(rowMeans(log(dhigh$mempdf))) # apply(-log(dhigh$mempdf), 1, mean)
dlow$favg <- exp(rowMeans(log(dlow$mempdf))) # apply(-log(dlow$mempdf), 1, mean)

dhlong <- dhigh %>%
  ungroup() %>%
  mutate(mempdf = as_tibble(mempdf)) %>%
  unnest(mempdf) %>%
  pivot_longer(c(enslin, dens, V1:V5, favg), names_to = "ensemble", values_to = "dens")

dllong <- dlow %>%
  ungroup() %>%
  mutate(mempdf = as_tibble(mempdf)) %>%
  unnest(mempdf) %>%
  pivot_longer(c(enslin, dens, V1:V5, favg), names_to = "ensemble", values_to = "dens")

p1 <- ggplot(test, aes(x = log_popularity, y = log_revenue)) +
  geom_point(alpha = 0.3) +
  geom_point(aes(y = response), data = dhigh, color = col_high) +
  geom_point(aes(y = response), data = dlow, color = col_low) +
  geom_hline(yintercept = dhigh$response[1], color = col_high, linetype = 5) +
  geom_hline(yintercept = dlow$response[1], color = col_low, linetype = 5) +
  theme_bw() + labs(x = expression(log[10]~popularity~score),
                    y = expression(log[10]~movie~revenue))

p2 <- ggplot(dhlong %>% filter(ensemble %in% c("dens", "enslin")), aes(x = log_revenue, y = dens, linetype = ensemble, linewidth = ensemble)) +
  geom_line(aes(color = "high")) +
  geom_line(data = dllong %>% filter(ensemble %in% c("dens", "enslin")), aes(color = "low")) +
  geom_line(data = dhlong %>% filter(ensemble %in% paste0("V", 1:5)), aes(color = "high", linetype = "members")) +
  geom_line(data = dllong %>% filter(ensemble %in% c("dens", "enslin")), aes(color = "low")) +
  geom_line(data = dllong %>% filter(ensemble %in% paste0("V", 1:5)), aes(color = "low", linetype = "members")) +
  scale_color_manual(values = c(col_high, col_low)) +
  theme_bw() +
  labs(y = expression(Conditional~density~of~Y~"|"~X==x),
       x = expression(log[10]~movie~revenue),
       color = "Observation", linetype = element_blank()) +
  geom_vline(xintercept = c(dlow$response[1], dhigh$response[1]),
             col = c(col_low, col_high), linetype = 5) +
  scale_linetype_manual(labels = c("dens" = "TRF", "enslin" = "LIN", "members" = "members"),
                                   # "V1" = "members"), # , "V2" = "members", "V3" = "members",
                                   # "V4" = "members", "V5" = "members"),
                        values = c(1, 2, 3)) + # , rep(3, 5))) +
  scale_linewidth_manual(values = c(mlwd, mlwd, rep(0.5, 5))) +
  guides(linewidth = "none", linetype = "none")

p4 <- ggplot(dhlong %>% filter(ensemble %in% c("dens", "enslin", "favg")),
             aes(x = log_revenue, y = -log(dens), linetype = ensemble, linewidth = ensemble)) +
  geom_line(aes(color = "high")) +
  geom_line(data = dllong %>% filter(ensemble %in% c("dens", "enslin", "favg")), aes(color = "low")) +
  geom_line(data = dhlong %>% filter(ensemble %in% paste0("V", 1:5)), aes(color = "high", linetype = "members")) +
  geom_line(data = dllong %>% filter(ensemble %in% c("dens", "enslin", "favg")), aes(color = "low")) +
  geom_line(data = dllong %>% filter(ensemble %in% paste0("V", 1:5)), aes(color = "low", linetype = "members")) +
  scale_color_manual(values = c(col_high, col_low)) +
  theme_bw() +
  labs(y = expression(Negative~log~density~of~Y~"|"~X==x),
       x = expression(log[10]~movie~revenue),
       color = "Observation", linetype = element_blank()) +
  geom_vline(xintercept = c(dlow$response[1], dhigh$response[1]),
             col = c(col_low, col_high), linetype = 5) +
  scale_linetype_manual(labels = c("dens" = "TRF", "enslin" = "LIN", "favg" = "AVG", "members" = "members"),
                                   # "V1" = "members"), # , "V2" = "members", "V3" = "members",
                                   # "V4" = "members", "V5" = "members"),
                        values = c(1, 2, 4, 3)) + # , rep(3, 5))) +
  scale_linewidth_manual(values = c(mlwd, mlwd, mlwd, rep(0.5, 5))) +
  guides(linewidth = "none") + theme(legend.position = "top")

cl <- apply(
  tr <- do.call("cbind", predict(ens, newdata = test, type = "cdf")), 1, mean
)

p3 <- data.frame(cll = qlogis(cl), trr = qlogis(tr)) %>%
  pivot_longer(trr.1:trr.5, names_to = "member", values_to = "trafo") %>%
  ggplot(aes(x = cll, y = trafo, color = member)) +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(alpha = 0.3) +
  theme_bw() +
  scale_color_brewer(palette = "Dark2") +
  labs(y = expression(Transformation~h[m](y~"|"~x)),
       x = expression(Classical~ensemble~{F^{-1}}[Z]({bar(F)^c}[M](y~"|"~x))))

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend.grob = get_legend(p4))
ggsave(file.path(bpath, "figure4.pdf"), height = 6, width = 7)
