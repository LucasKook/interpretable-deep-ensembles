# Functions
# Andrea Goetschi
# April 2022


# Functions to create plots -----------------------------------------------

## Performance plots (test error)

pl_met <- function(spl_met,
                   indiv_met = NULL,
                   ci,
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
                                "cimrsbin" = expression(paste(CI[B], "-Binary")),
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
                   hjust = 0.5,
                   ebarwidth = 0.4,
                   refebarwidth = 0.1,
                   t.size = NULL) {
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
  # same levels
  if (!is.factor(spl_met$mod)) spl_met$mod <- factor(spl_met$mod)
  ci <- ci %>% filter(mod %in% unique(spl_met$mod)) %>%
    mutate(mod = factor(mod, levels = levels(spl_met$mod)))

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
    (if(!is.null(indiv_met)) geom_beeswarm(inherit.aes = FALSE,
                                           data = indiv_met %>% filter(metric %in% metrics),
                                           aes(x = mod, y = val), 
                                           alpha = 0.2, size = 0.3, dodge.width = 1)) +
    stat_summary(fun = mean, geom = "point", position = position_dodge2(width = 0.75),
                 shape = "diamond",
                 size = 1.8) +
    geom_errorbar(data = ci %>% {if (!rel) filter(., metric %in% metrics) %>%
                                           mutate(metric = factor(metric, levels = metrics))
                                 else filter(., metric %in% paste0("d", metrics)) %>%
                                      mutate(metric = factor(sub('.', '', metric), levels = metrics))}, 
                  aes(x = mod, ymin = lwr, ymax = upr,
                      group = interaction(mod, method),
                      color = method), 
                  position = position_dodge(width = 0.75),
                  size = 0.3, width = ebarwidth,
                  inherit.aes = FALSE) +
    (if (!is.null(ref)) stat_summary(inherit.aes = FALSE,
                                     fun = mean,
                                     geom = "point",
                                     shape = "diamond",
                                     data = spl_met %>% filter(mod %in% ref, metric %in% metrics),
                                     aes(x = mod, y = val),
                                     color = col_ref, size = 1.8)) + 
    labs(x = xlab, y = if (!rel) ylab else parse(text = ylab_rel)) +
    scale_color_manual(name = legend_title,
                       labels = if (!is.null(ref))
                                  ens_labs[names(ens_labs) %in% levels(spl_met$method)[levels(spl_met$method) != unique(spl_met[spl_met$mod %in% ref, "method"])]]
                                else ens_labs_noref,
                                  values = col_ens,
                       guide = guide_legend(byrow = FALSE,
                                            reverse = TRUE),
                       breaks = if (!is.null(ref))
                                  levels(spl_met$method)[levels(spl_met$method) != unique(spl_met[spl_met$mod %in% ref, "method"])]
                                else levels(spl_met$method)
    ) +
    scale_shape_discrete(guide = "none") +
    scale_x_discrete(labels = mod_labs, limits = levels(spl_met$mod)) +
    theme_bw() +
    (if (!is.null(t.size)) theme(text = element_text(size = t.size))) +
    theme(axis.text.x = element_text(angle = angle, vjust = vjust, hjust = hjust),
          panel.grid.major.y = element_blank(),
          legend.position = legend_pos,
          strip.background = element_rect(fill = "white")) +
    (if (!ticks) theme(axis.text.y = element_blank(), 
                       axis.ticks.y = element_blank())) +
    (if (!legend) theme(legend.position = "none")) +
    coord_flip()
}


## Plot log odds-ratios

pl_or <- function(indiv,
                  members = TRUE,
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
                               "x_10" = expression(x[10]),
                               "mrs_before4" = "mRS at BL 4", 
                               "mrs_before2" = "mRS at BL 2",
                               "nihss_baseline" = "NIHSS at BL", 
                               "mrs_before3" = "mRS at BL 3",
                               "mrs_before1" = "mRS at BL 1", 
                               "rf_hypertoniay" = "hypertension",
                               "stroke_beforey" = "prior stroke", 
                               "rf_smokery" = "smoking",
                               "rf_chdy" = "CHD", 
                               "rf_atrial_fibrillationy" = "atrial fibrillation",        
                               "age" = "age", 
                               "tia_beforey" = "prior TIA", 
                               "rf_diabetesy" = "diabetes",
                               "sexm" = "male",
                               "rf_hypercholesterolemiay" = "hypercholesterolemia"),
                  mod_labs = c("silscs" = expression(paste("SI-", CS[B], "-", LS[x])),
                               "cilsmrsbl" = expression(paste(CI[B], "-", LS[mRS])),
                               "cils" = expression(paste(CI[B], "-", LS[x])),
                               "sils" = expression(paste("SI-", LS[x]))),
                  mod_cols = c("sils" = "grey40",
                               "cils" =  qualitative_hcl(n = 3, l = 40)[1],
                               "cilsmrsbl" = qualitative_hcl(n = 3, l = 40)[2],
                               "silscs" =  qualitative_hcl(n = 3, l = 40)[3]),
                  width = 0.6,
                  ylim = NULL,
                  e = FALSE,
                  weighted = TRUE,
                  avg_across_spl = TRUE,
                  refline = FALSE,
                  y.text = TRUE,
                  t.size = NULL,
                  lbetvar = FALSE,
                  pooled_only = FALSE,
                  ci_sils = NULL,
                  order_vars = TRUE) {
  
  uw <- indiv %>% group_by(var, mod, spl) %>% summarise(avg = mean(lor),
                                                        lwr = bs(lor, weights = rep(1, length(lor)))[1],
                                                        upr = bs(lor, weights = rep(1, length(lor)))[2])
  uw <- uw %>% mutate(weights = "equal")
  # 1000 values per model, var, spl
  uw_pooled <- indiv %>% filter(mod != "sils") %>% group_by(var, mod, spl) %>% summarise(bs = bs(lor, ret_all = T)) %>%
    group_by(var, mod) %>% summarise(lwr = quantile(bs, 0.025, na.rm = T),
                                     upr = quantile(bs, 0.975, na.rm = T))
  
  if (is.null(ci_sils)) {
    uw_pooled_sils <- indiv %>% filter(mod == "sils") %>% group_by(var, mod) %>% summarise(bs = bs(lor, ret_all = T)) %>%
      group_by(var, mod) %>% summarise(lwr = quantile(bs, 0.025, na.rm = T),
                                       upr = quantile(bs, 0.975, na.rm = T))
    uw_pooled <- bind_rows(uw_pooled, uw_pooled_sils)
  } else {
    uw_pooled <- bind_rows(uw_pooled, ci_sils)
  }
  
  uw_pooled <- uw_pooled %>% mutate(weights = "equal")
  uw_pooled_avg <- uw %>% group_by(var, mod) %>% summarise(avg = mean(avg))
  uw_pooled <- left_join(uw_pooled, uw_pooled_avg, by = c("var", "mod"))
  
  if (!weighted) {
    spl_avg <- uw
    pooled <- uw_pooled
  } else {
    w <- indiv %>% group_by(var, mod, spl) %>% summarise(avg = weighted.mean(lor, w = w),
                                                         lwr = bs(lor, weights = w)[1],
                                                         upr = bs(lor, weights = w)[2])
    w <- w %>% mutate(weights = "tuned")
    
    w_pooled <- indiv %>% filter(mod != "sils") %>% group_by(var, mod, spl) %>% summarise(bs = bs(lor, weights = w, ret_all = T)) %>%
      group_by(var, mod) %>% summarise(lwr = quantile(bs, 0.025, na.rm = T),
                                       upr = quantile(bs, 0.975, na.rm = T))
    if (is.null(ci_sils)) {
      w_pooled <- bind_rows(w_pooled, uw_pooled_sils) # SI_LS: no weights
    } else {
      w_pooled <- bind_rows(w_pooled, ci_sils)
    }
    w_pooled <- w_pooled %>% mutate(weights = "tuned")
    w_pooled_avg <- w %>% group_by(var, mod) %>% summarise(avg = mean(avg))
    w_pooled <- left_join(w_pooled, w_pooled_avg, by = c("var", "mod"))
    
    spl_avg <- bind_rows(uw, w)
    pooled <- bind_rows(uw_pooled, w_pooled)
  }
  if (e) {
    spl_avg$avg <- exp(spl_avg$avg)
    indiv$lor <- exp(indiv$lor)
    pooled$lwr <- exp(pooled$lwr)
    pooled$upr <- exp(pooled$upr)
    pooled$avg <- exp(pooled$avg)
  }
  
  # order levels of vars
  if (order_vars) {
    vord <- indiv %>% group_by(var) %>% summarise(avg = mean(lor)) %>% arrange(avg)
    vord <- unique(vord$var)
    indiv$var <- factor(indiv$var, levels = vord)
    spl_avg$var <- factor(spl_avg$var, levels = vord)
    pooled$var <- factor(pooled$var, levels = vord)
  }
  
  # order models 
  spl_avg$mod <- factor(spl_avg$mod, levels = levels(indiv$mod))
  pooled$mod <- factor(pooled$mod, levels = levels(indiv$mod))
  
  if (!pooled_only) {
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
      
      (if (avg_across_spl) geom_point(data = pooled, 
                                      aes(x = var, y = avg,
                                          group = interaction(var, mod), color = mod),
                                      # position = position_dodge(width = 0.7),
                                      position = position_dodgenudge(x = -0.48, y = 0, width = 0.1),
                                      shape = "diamond",
                                      size = 2.5, alpha = 0.4, inherit.aes = FALSE)) +
      (if (avg_across_spl) geom_errorbar(data = pooled,
                                         aes(x = var, ymin = lwr, ymax = upr,
                                             group = interaction(var, mod), color = mod),
                                         # position = position_dodge(width = 0.7),
                                         position = position_dodgenudge(x = -0.48, y = 0, width = 0.5),
                                         width = 0.1,
                                         size = 0.8, alpha = 0.4, inherit.aes = FALSE)) +
      # geom_beeswarm(dodge.width = 0.75, size = 0.5, alpha = 0.7) +
      (if (members) geom_beeswarm(data = indiv %>% filter(mod != "sils"),
                                  aes(x = var, y = lor,
                                      group = interaction(var, mod, spl)),
                                  dodge.width = 0.75, size = 0.3, alpha = 0.2, #dodge.width = 0.75
                                  inherit.aes = FALSE)) +
      (if (refline) geom_hline(yintercept = if (!e) 0 else 1, alpha = 0.2)) +
      # labs(x = "", y = parse(text = if(!e) "lower~risk%<-%hat(beta)%->%higher~risk"
      #                        else "lower~risk%<-%exp(hat(beta))%->%higher~risk")) +
      labs(x = "", y = parse(text = if(!e) "hat(beta)"
                             else "exp(hat(beta))")) +
      (if (!is.null(ylim)) lims(y = ylim)) +
      scale_color_manual(name = "Models",
                         labels = mod_labs,
                         values = mod_cols,
                         breaks = levels(indiv$mod)) +
      scale_x_discrete(labels = var_labs) +
      coord_flip() +
      theme_bw() + 
      (if (lbetvar) geom_vline(xintercept = seq(1.5, length(unique(indiv$var)), by = 1), alpha = 0.05)) +
      (if (!y.text) theme(axis.ticks.y = element_blank(),
                          axis.text.y = element_blank())) +
      (if (!is.null(t.size)) theme(text = element_text(size = t.size))) +
      theme(panel.grid.major.y = element_blank(),
            strip.background = element_rect(fill = "white")) 
  } else {
    ggplot(spl_avg, aes(x = var, y = avg, group = interaction(var, mod))) +
      (if(weighted) facet_grid(rows = vars(weights), labeller = labeller(.rows = label_both))) +
      (if (avg_across_spl) geom_point(data = pooled, 
                                      aes(x = var, y = avg,
                                          group = interaction(var, mod), color = mod),
                                      # position = position_dodge(width = 0.7),
                                      position = position_dodge(width = 0.75),
                                      shape = "diamond",
                                      size = 2, inherit.aes = FALSE)) +
      (if (avg_across_spl) geom_errorbar(data = pooled,
                                         aes(x = var, ymin = lwr, ymax = upr,
                                             group = interaction(var, mod), color = mod),
                                         position = position_dodge(width = 0.75), width = width,
                                         size = 0.5, inherit.aes = FALSE)) +
      (if (refline) geom_hline(yintercept = if (!e) 0 else 1, alpha = 0.2)) +
      # labs(x = "", y = parse(text = if(!e) "lower~risk%<-%hat(beta)%->%higher~risk"
      #                        else "lower~risk%<-%exp(hat(beta))%->%higher~risk")) +
      labs(x = "", y = parse(text = if(!e) "hat(beta)"
                             else "exp(hat(beta))")) +
      (if (!is.null(ylim)) lims(y = ylim)) +
      scale_color_manual(name = "Models",
                         labels = mod_labs,
                         values = mod_cols,
                         breaks = levels(indiv$mod)) +
      scale_x_discrete(labels = var_labs) +
      coord_flip() +
      theme_bw() + 
      (if (lbetvar) geom_vline(xintercept = seq(1.5, length(unique(indiv$var)), by = 1), alpha = 0.05)) +
      (if (!y.text) theme(axis.ticks.y = element_blank(),
                          axis.text.y = element_blank())) +
      (if (!is.null(t.size)) theme(text = element_text(size = t.size))) +
      theme(panel.grid.major.y = element_blank(),
            strip.background = element_rect(fill = "white")) 
  }
}


## Calibration plots

pl_cal <- function(avg, avg_ref = NULL, spl = NULL, spl_ref = NULL, indiv = NULL,
                   facet_labs = c("sils" = "SI*'-'*LS[x]",
                                  "si" = "SI",
                                  "sics" = "SI*'-'*CS[B]",
                                  "cils" = "CI[B]*'-'*LS[x]",
                                  "cilsmrsbl" = "CI[B]*'-'*LS[mRS]",
                                  "ci" = "CI[B]",
                                  "cimrsbinary" = "CI[B]*'-'*Binary",
                                  "silscs" = "SI*'-'*CS[B]*'-'*LS[x]"),
                   weighted = TRUE,
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
                   col_ref = "grey40",
                   psize = 0.7, lsize = 0.25, ebarsize = 0.25,
                   strokep = FALSE,
                   nrow = 2, ncol = 4,
                   t.size = NULL) {
  
  ggplot(avg, aes(x = midpoint, y = prop, color = method)) +
    (if (!strokep) {
    (if (weighted) { facet_grid(weights ~ mod, 
                              labeller = labeller(mod = as_labeller(facet_labs, label_parsed),
                                   .rows = label_both))
     } else {
       facet_grid(~ mod, labeller = labeller(mod = as_labeller(facet_labs, label_parsed)))
       })} else {
         facet_wrap(~ mod, labeller = labeller(mod = as_labeller(facet_labs, label_parsed)), 
                    ncol = ncol, nrow = nrow)
       })+
    geom_point(size = psize) +
    geom_line(size = lsize) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), size = ebarsize, width = 0.06) +
    (if(!is.null(avg_ref)) geom_point(data = avg_ref, aes(x = midpoint, y = prop), 
                                      color = col_ref, size = psize, inherit.aes = FALSE)) +
    (if(!is.null(avg_ref)) geom_line(data = avg_ref, aes(x = midpoint, y = prop), 
                                     color = col_ref, size = lsize, inherit.aes = FALSE)) +
    (if(!is.null(avg_ref)) geom_errorbar(data = avg_ref, aes(x = midpoint, ymin = lwr, ymax = upr), 
                                         color = col_ref, size = ebarsize, inherit.aes = FALSE,
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
    (if(!is.null(spl)) geom_point(data = spl, aes(x = midpoint, y = prop, 
                                  color = method), size = 0.7, alpha = 0.3, inherit.aes = FALSE)) +
    (if(!is.null(spl)) geom_line(data = spl, aes(x = midpoint, y = prop,
                                                 color = method,
                                                 group = interaction(spl, method)),
                                 alpha = 0.3, size = 0.25, inherit.aes = FALSE)) +
    (if(!is.null(spl)) geom_errorbar(data = spl, aes(x = midpoint, ymin = lwr, ymax = upr, 
                                                     color = method,
                                                     group = interaction(spl, method),
                                                     width = 0.06),
                                     alpha = 0.3, size = 0.1, inherit.aes = FALSE)) +
    (if(!is.null(spl_ref)) geom_point(data = spl_ref, aes(x = midpoint, y = prop, group = spl, 
                                      color = col_ref), size = 0.7, alpha = 0.3, inherit.aes = FALSE)) +
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
    (if (!is.null(t.size)) theme(text = element_text(size = t.size))) +
    (if (!legend) theme(legend.position = "none")) +
    theme(panel.grid.major.y = element_blank(),
          strip.background = element_rect(fill = "white"))
}



# Other functions ---------------------------------------------------------


# Create list of CDFs and true Y (all splits) -----------------------------

load_lys_cdf_all <- function(all, m, K, 
                             l = c("nll", "rps"),
                             t = c("test", "val", "train")) {
  l <- match.arg(l)
  t <- match.arg(t)
  lys_spl <- all %>% filter(mod == m, loss == l, type == t) %>%
    group_split(spl)
  ret <- lapply(lys_spl, function(x) {
    tmp <- x %>% group_split(mem)
    lapply(tmp, function(z) select(z, all_of(1:K)))
  })
  return(ret)
}

load_y_true_all <- function(all, K, t = c("test", "val", "train")) {
  t <- match.arg(t)
  tmp <- all %>% filter(type == t) %>% group_split(spl)
  ret <- lapply(tmp, function(x) select(x, all_of(1:K)) %>% as.matrix)
  return(ret)
}

# Bootstrap for log OR ----------------------------------------------------

bs <- function(vals, B = 1000, weights = rep(1, length(vals)), ret_all = FALSE) {
  n <- length(vals)
  avgs <- numeric(B)
  for (b in seq_len(B)) {
    idx <- sample(n, n, replace = TRUE)
    avgs[b] <- weighted.mean(vals[idx], weights[idx])
  }
  if (!ret_all) {
    return(quantile(avgs, c(0.025, 0.975), na.rm = TRUE)) # otherwise error (all w = 0 --> avgs = NaN)
  } else {
    return(avgs)
  }
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
                    custom_cuts_fun = function(x) quantile(x, c(0.5, 0.9, 0.99, 0.999))
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
                       custom_cuts_fun = function(x) quantile(x, c(0.5, 0.9, 0.99, 0.999))
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
                         custom_cuts_fun = function(x) quantile(x, c(0.5, 0.9, 0.99, 0.999))
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
                         custom_cuts_fun = function(x) quantile(x, c(0.5, 0.9, 0.99, 0.999))
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
                      custom_cuts_fun = function(x) quantile(x, c(0.5, 0.9, 0.99, 0.999))
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

extract_w <- function(w_all, meth = c("linear", "log-linear", "trafo"), ens = 5) {
  meth <- match.arg(meth)
  w_all %>% 
    filter(method  == meth) %>% 
    select(1:all_of(ens)) %>%
    t() %>% as.data.frame() %>% as.list()
}
