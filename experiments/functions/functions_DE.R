# Functions
# Andrea Goetschi
# April 2022


# Functions to create plots -----------------------------------------------

## Performance plots (test error)

pl_met <- function(spl_met,
                   indiv_met = NULL,
                   ref = NULL, 
                   rel = FALSE,
                   weighted = TRUE,
                   metrics = c("binnll", "brier", "eauc", "ebinacc"),
                   facet_labs = c("binnll" = "Binary NLL",
                                  "brier" = "Brier score",
                                  "eauc" = "1 - AUC",
                                  "ebinacc" = "1 - ACC",
                                  "eacc" = "1 - ACC",
                                  "nll" = "NLL",
                                  "rps" = "RPS",
                                  "eqwk" = "1 - QWK",
                                  "cint" = "CITL",
                                  "cslope" = "C slope"),
                   ens_labs = c("avg" = "LIN-Avg",
                                "avgll" = "LOG-Avg",
                                "avgtrf" = "TRF-Avg",
                                "linear" = "LIN-Ens",
                                "log-linear" = "LOG-Ens",
                                "trafo" = "TRF-Ens",
                                "ref" = "XX"),
                   ens_labs_noref = c("avg" = "LIN-Avg",
                                      "avgll" = "LOG-Avg",
                                      "avgtrf" = "TRF-Avg",
                                      "linear" = "LIN-Ens",
                                      "log-linear" = "LOG-Ens",
                                      "trafo" = "TRF-Ens"),
                   mod_labs = c("silscs" = expression(paste("SI-", CS[B], "-", LS[x])),
                                "cils" = expression(paste(CI[B], "-", LS[x])),
                                "cilsmrsbl" = expression(paste(CI[B], "-", LS[mRS])),
                                "sics" = expression(paste("SI-", CS[B])),
                                "ci" = expression(paste(CI[B])),
                                "sils" = expression(paste("SI-", LS[x])),
                                "cimrsbin" = "MCC-Binary",
                                "si" = "SI"),
                   legend = TRUE,
                   legend_title = "Method",
                   legend_pos = "right",
                   xlab = expression(lower%<-%complexity%->%higher),
                   ylab = expression(paste(better%<-%test, " ", error%->%worse)),
                   ylab_rel = "SI*'-'*LS[x]~worse%<-%difference~'in'~test~error%->%SI*'-'*LS[x]~better",
                   ticks = TRUE,
                   col_ens = c("avg" = qualitative_hcl(n = 3, l = 80)[1],
                               "avgll" = qualitative_hcl(n = 3, l = 80)[2],
                               "avgtrf" = qualitative_hcl(n = 3, l = 80)[3],
                               "linear" = qualitative_hcl(n = 3, l = 40)[1],
                               "log-linear" = qualitative_hcl(n = 3, l = 40)[2],
                               "trafo" = qualitative_hcl(n = 3, l = 40)[3]),
                   col_ref = "grey40",
                   nrow = 1, ncol = 4,
                   angle = 0,
                   vjust = 0.5,
                   hjust = 0.5) {
  if (rel) {
    sval <- spl_met %>%
      {if (weighted) filter(., mod == "sils", weights == "equal") else filter(., mod == "sils")} %>%
      dplyr::select(val, spl, metric)
    
    spl_met <- left_join(spl_met, sval, by = c("metric", "spl"))
    
    spl_met <- spl_met %>% 
      {if (weighted) group_by(., mod, method, weights) else group_by(., mod, method)} %>%
      mutate(val = val.x - val.y) %>%
      filter(mod != "sils")
    
    spl_met$mod <- factor(spl_met$mod, 
                          levels = levels(spl_met$mod)[!(levels(spl_met$mod) %in% "sils")])
  }
  ggplot(data = spl_met %>% filter(metric %in% metrics, !(mod %in% ref)), 
         aes(x = mod, y = val, group = interaction(mod, method),
             color = method, shape = spl)) +
    (if (weighted) facet_grid(weights ~ metric, scales = "free_x", 
                              labeller = labeller(.cols = facet_labs, .rows = label_both)) 
     else facet_grid(~ metric, scales = "free_x", 
                     labeller = labeller(.cols = facet_labs, .rows = label_both))) +
    geom_beeswarm(dodge.width = 0.75, size = 0.7, alpha = 0.55) +
    (if (rel) geom_hline(yintercept = 0, alpha = 0.2)) +
    (if ("cint" %in% metrics) geom_hline(data = data.frame(metric = "cint",
                                                           val = 0),
                                         aes(yintercept = val), alpha = 0.2)) +
    (if ("cslope" %in% metrics) geom_hline(data = data.frame(metric = "cslope",
                                                             val = 1),
                                           aes(yintercept = val), alpha = 0.2)) +
    (if (!is.null(ref)) geom_beeswarm(inherit.aes = FALSE,
                                      dodge.width = 0.5,
                                      size = 0.7, col = col_ref, alpha = 0.55,
                                      data = spl_met %>% filter(mod %in% ref, metric %in% metrics),
                                      aes(x = mod, y = val, shape = spl))) +
    # (if(!is.null(indiv_met)) geom_rug(inherit.aes = FALSE,
    #                                   data = indiv_met %>% filter(metric %in% metrics),
    #                                   aes(y = val), alpha = 0.4)) +
    (if(!is.null(indiv_met)) geom_beeswarm(inherit.aes = FALSE,
                                           data = indiv_met %>% filter(metric %in% metrics),
                                           aes(x = mod, y = val), 
                                           alpha = 0.2, size = 0.3, dodge.width = 1)) +
    # stat_summary(fun.data = "mean_cl_boot", position = position_dodge2(width = 0.75), size = 0.3) +
    stat_summary(fun = mean, geom = "point", position = position_dodge2(width = 0.75),
                 shape = "diamond",
                 size = 1.8) +
    stat_summary(fun.data = "mean_cl_boot", position = position_dodge(width = 0.75), 
                 size = 0.3, geom = "errorbar", width = 0.4) +
    (if (!is.null(ref)) stat_summary(inherit.aes = FALSE,
                                     fun = mean,
                                     geom = "point",
                                     shape = "diamond",
                                     data = spl_met %>% filter(mod %in% ref, metric %in% metrics),
                                     aes(x = mod, y = val),
                                     color = col_ref, size = 1.8)) + 
    (if (!is.null(ref)) stat_summary(inherit.aes = FALSE,
                                     fun.data = "mean_cl_boot",
                                     geom = "errorbar", width = 0.2,
                                     data = spl_met %>% filter(mod %in% ref, metric %in% metrics),
                                     aes(x = mod, y = val),
                                     color = col_ref, size = 0.3)) + 
    labs(x = xlab, y = if (!rel) ylab else parse(text = ylab_rel)) +
    (if (legend) scale_color_manual(name = legend_title,
                                    labels = if (!is.null(ref))
                                      ens_labs[names(ens_labs) %in% levels(spl_met$method)[levels(spl_met$method) != unique(spl_met[spl_met$mod %in% ref, "method"])]]
                                    else ens_labs_noref,
                                    values = col_ens,
                                    guide = guide_legend(byrow = FALSE,
                                                         reverse = TRUE),
                                    breaks = if (!is.null(ref))
                                      levels(spl_met$method)[levels(spl_met$method) != unique(spl_met[spl_met$mod %in% ref, "method"])]
                                    else levels(spl_met$method)
    )) +
    scale_shape_discrete(guide = "none") +
    scale_x_discrete(labels = mod_labs, limits = levels(spl_met$mod)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust = hjust),
          panel.grid.major.y = element_blank(), # aspect.ratio = 1,
          legend.position = legend_pos,
          strip.background = element_rect(fill = "white")) +
    (if (!ticks) theme(axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank())) +
    (if (!legend) theme(legend.position = "none")) +
    coord_flip()
}


## Plot log odds-ratios

pl_or <- function(indiv,
                  weights,
                  var_labs = c(gender = "female",
                               "x_1" = expression(x[1]),
                               "x_2" = expression(x[2]),
                               "x_3" = expression(x[3]),
                               "x_4" = expression(x[4]),
                               "x_5" = expression(x[5]),
                               "x_6" = expression(x[6]),
                               "x_7" = expression(x[7]),
                               "x_8" = expression(x[8]),
                               "x_9" = expression(x[9]),
                               "x_10" = expression(x[10])),
                  mod_labs = c("silscs" = expression(paste("SI-", CS[B], "-", LS[x])),
                               "cilsmrsbl" = expression(paste(CI[B], "-", LS[mRS])),
                               "cils" = expression(paste(CI[B], "-", LS[x])),
                               "sils" = expression(paste("SI-", LS[x]))),
                  mod_cols = c("sils" = "grey40",
                               "cils" =  qualitative_hcl(n = 2, l = 40)[1],
                               "silscs" =  qualitative_hcl(n = 2, l = 40)[2]),
                  width = 0.6,
                  ylim = NULL,
                  e = FALSE,
                  weighted = TRUE,
                  avg_across_spl = TRUE) {
  
  uw <- indiv %>% group_by(var, mod, spl) %>% summarise(avg = mean(lor),
                                                        lwr = bs(lor, weights = rep(1, length(lor)))[1],
                                                        upr = bs(lor, weights = rep(1, length(lor)))[2])
  uw <- uw %>% mutate(weights = "equal")
  if (!weighted) {
    spl_avg <- uw
  } else {
    w <- indiv %>% group_by(var, mod, spl) %>% summarise(avg = weighted.mean(lor, w = w),
                                                         lwr = bs(lor, weights = w)[1],
                                                         upr = bs(lor, weights = w)[2])
    w <- w %>% mutate(weights = "tuned")
    spl_avg <- bind_rows(uw, w)
  }
  if (e) {
    spl_avg$avg <- exp(spl_avg$avg)
    indiv$lor <- exp(indiv$lor)
  }
  if (avg_across_spl) {
    avg_acr_spl <- spl_avg %>% group_by(var, mod, weights) %>% summarise(avgspl = mean(avg),
                                                                            lwr = bs(avg, weights = rep(1, length(avg)))[1],
                                                                            upr = bs(avg, weights = rep(1, length(avg)))[2])
  }
  
  ggplot(spl_avg, aes(x = var, y = avg, group = interaction(var, mod, spl))) +
    (if(weighted) facet_grid(rows = vars(weights), labeller = labeller(.rows = label_both))) +
    geom_point(aes(x = var, y = avg,
                   group = interaction(var, mod, spl), color = mod),
               shape = "diamond", position = position_dodge(width = 0.75),
               size = 2) +
    geom_errorbar(aes(ymin = lwr, ymax = upr,
                      group = interaction(var, mod, spl), color = mod),
                  position = position_dodge(width = 0.75), width = width,
                  size = 0.5) +
    
    (if (avg_across_spl) geom_point(data = avg_acr_spl, 
                                    aes(x = var, y = avgspl,
                                        group = interaction(var, mod), color = mod),
                                    position = position_dodge(width = 0.3),
                                    shape = "diamond",
                                    size = 3, alpha = 0.2, inherit.aes = FALSE)) +
    (if (avg_across_spl) geom_errorbar(data = avg_acr_spl,
                                       aes(x = var, ymin = lwr, ymax = upr,
                                           group = interaction(var, mod), color = mod),
                                       position = position_dodge(width = 0.3),
                                       width = 0.2,
                                       size = 1, alpha = 0.2, inherit.aes = FALSE)) +
    
    # geom_beeswarm(dodge.width = 0.75, size = 0.5, alpha = 0.7) +
    geom_beeswarm(data = indiv,
                  aes(x = var, y = lor,
                      group = interaction(var, mod, spl)),
                  dodge.width = 0.75, size = 0.3, alpha = 0.2,
                  inherit.aes = FALSE) +
    geom_hline(yintercept = if (!e) 0 else 1, alpha = 0.2) +
    labs(x = "", y = parse(text = if(!e) "lower~risk%<-%hat(beta)%->%higher~risk"
                           else "lower~risk%<-%exp(hat(beta))%->%higher~risk")) +
    (if (!is.null(ylim)) lims(y = ylim)) +
    scale_color_manual(name = "Models",
                       labels = mod_labs,
                       values = mod_cols,
                       breaks = levels(indiv$mod)) +
    scale_x_discrete(labels = var_labs) +
    coord_flip() +
    theme_bw() + 
    theme(panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "white"))
}


## Calibration plots

pl_cal <- function(avg, avg_ref = NULL, spl = NULL, spl_ref = NULL, indiv = NULL,
                   facet_labs = c("sils" = "SI*'-'*LS[x]",
                                  "si" = "SI",
                                  "sics" = "SI*'-'*CS[B]",
                                  "cils" = "CI[B]*'-'*LS[x]",
                                  "cilsmrsbl" = "CI[B]*'-'*LS[mRS]",
                                  "ci" = "CI[B]",
                                  "cimrsbinary" = "MCC-Binary",
                                  "silscs" = "SI*'-'*CS[B]*'-'*LS[x]"),
                   legend = TRUE,
                   legend_title = "Method",
                   ens_labs = c("avg" = "LIN-Avg",
                                "avgll" = "LOG-Avg",
                                "avgtrf" = "TRF-Avg",
                                "linear" = "LIN-Ens",
                                "log-linear" = "LOG-Ens",
                                "trafo" = "TRF-Ens"),
                   col_ens = c("avg" = qualitative_hcl(n = 3, l = 80)[1],
                               "avgll" = qualitative_hcl(n = 3, l = 80)[2],
                               "avgtrf" = qualitative_hcl(n = 3, l = 80)[3],
                               "linear" = qualitative_hcl(n = 3, l = 40)[1],
                               "log-linear" = qualitative_hcl(n = 3, l = 40)[2],
                               "trafo" = qualitative_hcl(n = 3, l = 40)[3]),
                   col_ref = "grey40") {
  
  ggplot(avg, aes(x = midpoint, y = prop, color = method)) +
    facet_grid(weights ~ mod, 
               labeller = labeller(mod = as_labeller(facet_labs, label_parsed),
                                   .rows = label_both)) +
    geom_point(size = 0.7) +
    geom_line(size = 0.25) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), size = 0.25, width = 0.06) +
    (if(!is.null(avg_ref)) geom_point(data = avg_ref, aes(x = midpoint, y = prop), 
                                      color = col_ref, size = 0.7, inherit.aes = FALSE)) +
    (if(!is.null(avg_ref)) geom_line(data = avg_ref, aes(x = midpoint, y = prop), 
                                     color = col_ref, size = 0.25, inherit.aes = FALSE)) +
    (if(!is.null(avg_ref)) geom_errorbar(data = avg_ref, aes(x = midpoint, ymin = lwr, ymax = upr), 
                                         color = col_ref, size = 0.25, inherit.aes = FALSE,
                                         width = 0.06)) +
    # geom_point(data = spl, aes(x = midpoint, y = prop,
    #                                color = method,
    #                                group = interaction(spl, method)),
    #            alpha = 0.2, size = 0.25, inherit.aes = FALSE) +
    # geom_point(data = spl_ref, aes(x = midpoint, y = prop,
    #                            group = spl),
    #            alpha = 0.2, size = 0.25, inherit.aes = FALSE) +
    (if (!is.null(indiv)) geom_line(data = indiv, aes(x = midpoint, y = prop,
                                                      group = interaction(ens, spl)),
                                    alpha = 0.2, size = 0.25, inherit.aes = FALSE)) +
    (if (!is.null(indiv)) geom_errorbar(data = indiv, aes(x = midpoint, ymin = lwr, ymax = upr,
                                                          group = interaction(ens, spl),
                                                          width = 0.06),
                                        alpha = 0.2, size = 0.25, inherit.aes = FALSE)) +
    (if(!is.null(spl)) geom_line(data = spl, aes(x = midpoint, y = prop,
                                                 color = method,
                                                 group = interaction(spl, method)),
                                 alpha = 0.3, size = 0.25, inherit.aes = FALSE)) +
    (if(!is.null(spl)) geom_errorbar(data = spl, aes(x = midpoint, ymin = lwr, ymax = upr, 
                                                     color = method,
                                                     group = interaction(spl, method),
                                                     width = 0.06),
                                     alpha = 0.3, size = 0.1, inherit.aes = FALSE)) +
    (if(!is.null(spl_ref)) geom_line(data = spl_ref, aes(x = midpoint, y = prop,
                                                         group = spl),
                                     alpha = 0.3, size = 0.25, inherit.aes = FALSE)) +
    (if(!is.null(spl_ref)) geom_errorbar(data = spl_ref, aes(x = midpoint, ymin = lwr, ymax = upr, 
                                                             group = spl, width = 0.06),
                                         alpha = 0.3, size = 0.1, inherit.aes = FALSE)) +
    geom_abline(slope = 1, intercept = 0, alpha = 0.2) +
    labs(y = "observed proportion", x = "predicted probability") +
    lims(x = 0:1, y = 0:1) +
    scale_color_manual(name = legend_title,
                       labels = ens_labs,
                       values = col_ens,
                       # guide = guide_legend(byrow = FALSE,
                       #                      reverse = TRUE),
                       breaks = levels(avg$method)) +
    theme_bw() +
    (if (!legend) theme(legend.position = "none")) +
    theme(panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "white"))
}



# Other functions ---------------------------------------------------------

# Bootstrap ---------------------------------------------------------------

bs <- function(vals, B = 1000, weights = rep(1, length(vals))) {
  n <- length(vals)
  avgs <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample(n, n, replace = TRUE)
    avgs[b] <- weighted.mean(vals[idx], weights[idx])
  }
  return(quantile(avgs, c(0.05, 0.95), na.rm = TRUE)) # otherwise error (all w = 0 --> avgs = NaN)
} 

# Bind rows matching a specific pattern -----------------------------------

bindr <- function(pat1 = "met", pat2 = NULL) {
  obj <- ls(pattern = pat1, envir = .GlobalEnv)
  if (!is.null(pat2)) {
    obj <- grep(pattern = pat2, x = obj, value = TRUE)
  }
  ret <- bind_rows(mget(obj, envir = .GlobalEnv))
  return(ret)
}

# Specify order of n first levels -----------------------------------------

relev <- function(d, var, levs) {
  d[, var] <- factor(d[, var], 
                     levels = c(levs, 
                                unique(d[, var])[!(unique(d[, var]) %in% levs)]))
  return(d)
}

# Calibration -------------------------------------------------------------

# calibration for one cdf (1 ensemble or 1 ensemble member)
cal_ens <- function(cdf, y_true, cuts = 11, cumulative = TRUE, pool = TRUE, emp = FALSE,
                    custom_cuts_fun = function(x) c(0, quantile(x, c(0, 0.5, 0.9, 0.99, 0.999, 1))[2:5])
                    ) {
  cdf <- as.matrix(cdf)
  K <- ncol(cdf)
  pdf <- t(apply(cbind(0, cdf), 1, diff))
  lys_cal <- list()
  start <- 1
  if (cumulative) end <- K-1 else end <- K
  #if (K == 2) start <- end <- 2
  for (cl in start:end) {
    if (cumulative) {
      yt <- factor(apply(y_true, 1, function(x) sum(x[start:cl])), levels = c(0, 1)) #levels = c(1, 0)) P(Y<y_k)
      yp <- 1 - apply(pdf, 1, function(x) sum(x[start:cl]))
    } else {
      yt <- factor(y_true[, cl], levels = c(0, 1)) #levels = c(1, 0))
      yp <- 1 - pdf[, cl]
    }
    df <- data.frame(true = yt, pred = yp)
    if (emp) {
      cuts <- custom_cuts_fun(yp)
    }
    if (K == 2) {
      lys_cal <- list(calibration(true ~ pred, data = df, cuts = cuts)$data)
    } else {
      lys_cal[[cl]] <- calibration(true ~ pred, data = df, cuts = cuts)$data
    }
  }
  tmp <- lapply(lys_cal, function(x) {
    as.matrix(x[, c("Percent", "Lower", "Upper", "Count", "midpoint")])
  })
  avg <- apply(simplify2array(tmp), 1:2, mean)
  
  if (pool) {
    ret <- data.frame(bin = lys_cal[[1]]$bin,
                      prop = avg[, "Percent"]/100,
                      lwr = avg[, "Lower"]/100,
                      upr = avg[, "Upper"]/100,
                      cases = avg[, "Count"],
                      midpoint = avg[, "midpoint"]/100)
  } else {
    ret <- lys_cal
  }
  return(ret)
}

# calculate calibration for each ensemble member --> take mean (weighted or non-weighted)
cal_splavg <- function(lys_cdf, y_true, cumulative = TRUE,  cuts = 11,
                       weights = rep(1, length(lys_cdf)), emp = FALSE,
                       custom_cuts_fun = function(x) c(0, quantile(x, c(0, 0.5, 0.9, 0.99, 0.999, 1))[2:5])
                       ) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  lys_cal <- lapply(lys_cdf, function(x) {
    cal_ens(cdf = x, y_true = y_true, cumulative = cumulative, cuts = cuts, emp = emp,
            custom_cuts_fun = custom_cuts_fun)
  })
  tmp <- lapply(lys_cal, function(x) {
    as.matrix(x[, c("prop", "lwr", "upr", "cases", "midpoint")])
  })
  avg <- apply(simplify2array(tmp), 1:2, weighted.mean, w = weights)
  
  ret <- data.frame(bin = lys_cal[[1]]$bin,
                    prop = avg[, "prop"],
                    lwr = avg[, "lwr"],
                    upr = avg[, "upr"],
                    cases = avg[, "cases"],
                    midpoint = avg[, "midpoint"])
  return(ret)
}

# calibration for all ensembles (splits)
comb_ens_cal <- function(lys_cdf_ens, y_true_all, cumulative = TRUE, cuts = 11, emp = FALSE,
                         custom_cuts_fun = function(x) c(0, quantile(x, c(0, 0.5, 0.9, 0.99, 0.999, 1))[2:5])
                         ) {
  lys_cdf_ens <- lapply(lys_cdf_ens, as.matrix)
  lys <- mapply(function(ens, y_true, cuts, emp) {
    cal_ens(cdf = ens, y_true = y_true, cumulative = cumulative, cuts = cuts, emp = emp,
            custom_cuts_fun = custom_cuts_fun)
  }, ens = lys_cdf_ens, y_true = y_true_all, cuts = cuts, emp = emp,
  SIMPLIFY = FALSE)
  for (s in seq_len(length(lys))) {
    lys[[s]]$spl <- factor(s)
  }
  ret <- do.call("rbind", lys)
  return(ret)
}

# calibration of each member --> mean; for all splits
comb_avg_cal <- function(lys_cdf_all, y_true_all, cumulative = TRUE, cuts = 11,
                         weights = rep(list(rep(1, length(lys_cdf_all[[1]]))), length(lys_cdf_all)), emp = FALSE,
                         custom_cuts_fun = function(x) c(0, quantile(x, c(0, 0.5, 0.9, 0.99, 0.999))[2:5])
                         ) {
  lys <- mapply(function(lys_cdf, y_true, cuts, weights, emp) {
    cal_splavg(lys_cdf = lys_cdf, y_true = y_true, cumulative = cumulative, cuts = cuts, weights = weights, emp = emp,
               custom_cuts_fun = custom_cuts_fun)
  }, lys_cdf = lys_cdf_all, y_true = y_true_all, cuts = cuts, weights = weights, emp = emp,
  SIMPLIFY = FALSE)
  for (s in seq_len(length(lys))) {
    lys[[s]]$spl <- factor(s)
  }
  ret <- do.call("rbind", lys)
  return(ret)
}

# average calibration across all splits 
# separate for mod, weights, method
avg_across_spl <- function(lys_splitted) {
  lys_avg <- lapply(lys_splitted, function(x) {
    tmp <-  x %>% group_split(spl)
    tmp <- lapply(tmp, function(z) {
      as.matrix(z[, c("prop", "lwr", "upr", "cases", "midpoint")])
    })
    avg <- apply(simplify2array(tmp), 1:2, mean)
    #ci <- apply(simplify2array(tmp), 1:2, smean.cl.boot)
    data.frame(bin = if (length(unique(x$bin)) != length(avg[, "prop"])) NA else unique(x$bin),
               prop = avg[, "prop"],
               #lwr = ci["Lower", , "prop"],
               #upr = ci["Upper", , "prop"],
               lwr = avg[, "lwr"],
               upr = avg[, "upr"],
               cases = avg[, "cases"],
               midpoint = avg[, "midpoint"],
               method = unique(x$method),
               mod = unique(x$mod),
               weights = unique(x$weights))
  })
  ret <- do.call("rbind", lys_avg)
  return(ret)
}

cal_indiv <- function(lys_cdf, y_true, cumulative = TRUE, cuts = 11, emp = FALSE,
                      custom_cuts_fun = function(x) c(0, quantile(x, c(0, 0.5, 0.9, 0.99, 0.999, 1))[2:5])
                      ) {
  lys_cdf <- lapply(lys_cdf, as.matrix)
  tmp <- lapply(lys_cdf, cal_ens, y_true = y_true, cumulative = cumulative, cuts = cuts, emp = emp,
                custom_cuts_fun = custom_cuts_fun)
  for (e in seq_len(length(tmp))) {
    tmp[[e]]$ens <- factor(e)
  }
  ret <- do.call("rbind", tmp)
  return(ret)
}

# weighted ensemble for all splits
get_weighted_ens <- function(lys_cdf_all, lys_weigths, type = c("linear", "log-linear", "trafo")) {
  t <- match.arg(type)
  mapply(function(lys, w) {
    get_ensemble(lys_cdf = lys, weights = w, type = t)
  }, lys = lys_cdf_all, w = lys_weigths, SIMPLIFY = FALSE)
} 

extract_w <- function(w_all, meth = c("linear", "log-linear", "trafo"), ens = ensembles) {
  meth <- match.arg(meth)
  w_all %>% 
    filter(method  == meth) %>% 
    select(1:all_of(ensembles)) %>%
    t() %>% as.data.frame() %>% as.list()
}

recalibrate <- function(cdf, cdf_val, y_true_val) {
  cdf <- as.matrix(cdf)
  cdf_val <- as.matrix(cdf_val)
  pdf <- t(apply(cbind(0, cdf), 1, diff))[, 2]
  pdf_val <- t(apply(cbind(0, cdf_val), 1, diff))[, 2]
  logits <- qlogis(pdf)
  logits_val <- qlogis(pdf_val)
  yt_val <- factor(y_true_val[, 2], levels = c(0, 1))
  m <- glm(yt_val ~ logits_val, family = "binomial")
  preds_new <- predict(m, newdata = data.frame(logits_val = logits), type = "response")
  cdf_new <- data.frame(V1 = 1 - preds_new, V2 = 1)
  return(cdf_new)
}

save_recalibrated_cdf <- function(lys_cdf_all, lys_cdf_all_val, y_true_all_val, 
                                  fname, spl = 6, ens = 5, out_dir, ref = FALSE) {
  if (ref) {
    ens <- 1
  }
  for (s in seq_len(spl)) {
    y_true_val <- y_true_all_val[[s]]
    for (e in seq_len(ens)) {
      cdf <- lys_cdf_all[[s]][[e]]
      cdf_val <- lys_cdf_all_val[[s]][[e]]
      cdf_new <- recalibrate(cdf = cdf, cdf_val = cdf_val, y_true_val = y_true_val)
      if (ref) {
        write.csv(cdf_new, file = paste0(out_dir, fname, "_cdftest_recal_spl", s, ".csv"),
                  row.names = FALSE)
      } else {
        write.csv(cdf_new, file = paste0(out_dir, fname, "_cdftest_recal_spl", s, "_ens", e, ".csv"),
                  row.names = FALSE)
      }
    }
  }
}

read_recalibrated_cdf <- function(fname, in_dir, spl = 6, ens = 5, ref = FALSE) {
  f <- list.files(path = in_dir,
                  pattern = paste0(fname, "_cdftest_recal"))
  if (ref) ens <- 1
  ret <- rep(list(list()), spl)
  for (s in seq_len(spl)) {
    for (e in seq_len(ens)) {
      ret[[s]][[e]] <- read.csv(paste0(in_dir, f[(s - 1) * ens + e]))
    }
  }
  return(ret)
}
