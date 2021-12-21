make_trait_env_lms <- function(trait, env){
  # Need to turn the quoted trait and env names into symbolic names
  trait <- sym(trait)
  env <- sym(env)

  tXe_lm_upper  <- lm(eval(trait) ~ eval(env), data = p17_sitesum_env)
  return(tXe_lm_upper)
}
