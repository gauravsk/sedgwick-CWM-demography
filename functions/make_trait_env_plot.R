lm_eqn <- function(m){
  eq <- substitute(italic(p)~"="~pval~","~~italic(r)^2~"="~r2, 
                   list(r2 = format(summary(m)$r.squared, digits = 3),
                        pval = format(summary(m)$coefficients[2,4], digits = 3)))
  as.character(as.expression(eq));
}

make_trait_env_plot <- function(trait, env, lms_obj, xl, yl){
  # Need to turn the quoted trait and env names into symbolic names
  trait <- sym(trait)
  trait_sem <- sym(str_replace(trait,"_mean", replace = "_sem"))
  env <- sym(env)
  preddf <- p17_sitesum_env
  preddf <- preddf %>% mutate(upper = !!trait + 1.96*!!trait_sem,
                    lower = !!trait - 1.96*!!trait_sem)
  predictions <- predict(lms_obj, preddf , se.fit = T, interval = "confidence", level = 0.95)
  preddf <- cbind(preddf, fit = predictions$fit, se.fit = predictions$se.fit)
  
  pval <- summary(lms_obj)$coefficients[2,4]
  
  gg_upper <- 
    ggplot(preddf, aes(x = !!env, y = !!trait)) + 
    geom_point(size = 4, pch = 21, fill = alpha("black", .65), stroke = 1) + 
    geom_errorbar(aes(x = !!env, ymin = lower,
                      ymax = upper), size = .3, width = 0) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1,vjust = 1,
             label=lm_eqn(lms_obj),
             parse = TRUE) +
    xlab(xl) + ylab(yl) +
    #ecoevoapps::theme_apps() + 
    theme(axis.title = element_text(size = 12))
  if(pval < .05) {
    gg_upper <- gg_upper + 
      geom_ribbon(data = preddf, aes(x = !!env, ymin = fit.lwr, ymax = fit.upr),
                                     alpha = .2) +
      geom_line(data = preddf, aes(x = !!env, y = fit.fit), color = "darkred")
  }
  return(gg_upper)
}
