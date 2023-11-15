rm(list=ls())
library("bayesplot")
library("ggplot2")
library("rstan")
library(gridExtra)
# Extract fitted model parameters and predictions
setwd("D:/Project/CAMP/behav-models")
fit <- readRDS("./Results/fit_lognormal.rds")
stan_data <- readRDS("Data/Stan_ready_data.rds")
source("./utils.R")

# Trace plot
# stan_trace(fit,pars=c("R_mu_base","R_mu_delta")) #normal model
stan_trace(fit,pars="R_mu") #lognormal model

# Check for divergences
num_divergent <- sum(get_sampler_params(fit, inc_warmup = FALSE)[[1]][, "divergent__"])

# Posterior distribution
# stan_hist(fit,pars=c("R_mu_base","R_mu_delta")) #normal model
stan_hist(fit,pars="R_mu") #lognormal model

# autocorrelation and convergence
stan_ac(fit, pars="R_mu")
# stan_diag(fit,information="sample",chains=1)  # not work
stan_diag(fit,information="divergence",chain=1)
stan_diag(fit,information="treedepth",chain=1)
stan_diag(fit,information="stepsize",chain=1)
# stan_rhat(fit,pars="R_mu")
# stan_ess(fit,pars="R_mu")
# stan_mcse(fit, pars="R_mu")
#
# Extract Posterior samples
post_samples <- rstan::extract(fit,
                          pars = "R_mu",permuted = FALSE)
rm(fit); gc()

# organize posterior samples into [n_samples, n_chains] matrix
samples_matrix <- organize_samples(post_samples)
# compute Rhat and ESS
compute_rhat(samples_matrix$`R_mu[1,2]`)
compute_ess(samples_matrix$`R_mu[1,2]`)
# Rhat is the potential scale reduction factor, which is an estimate of the ratio of the between-chain variance
# to the within-chain variance. It is computed by running multiple chains and comparing the variances between chains
# to the variance within chains. If the chains have not mixed well, the between-chain variance will be larger than
# the within-chain variance (because all the chains are targeting different stationary distributions), so the Rhat
# will be larger than 1.1. If Rhat is close to 1, then the chains have mixed well and the within-chain variance is
# an accurate estimate of the true variance.

# ESS is the effective sample size, which is the number of independent samples that would contain the same amount of
# information as the posterior samples. It is computed by running multiple chains and comparing the variance within
# chains to the variance between chains. If the chains have not mixed well, the between-chain variance will be larger
# than the within-chain variance (because all the chains are targeting different stationary distributions), so the
# ESS will be smaller than the number of iterations. If the chains have mixed well, the ESS will be close to the
# number of iterations.
# Rhat and ESS are computed for each parameter in the model. If Rhat is close to 1 and ESS is close to the number of
# iterations, then the chains have mixed well and the samples are a good representation of the posterior distribution.
# If Rhat is much larger than 1 and ESS is much smaller than the number of iterations, then the chains have not mixed
# well and the samples are not a good representation of the posterior distribution.

# plot Bayesian result with Frequentist result
set.seed(43202)
n_subj <- stan_data$N
subjs <- sample(1:n_subj, 3, replace = F)
# combine all chains
post_samples <- rstan::extract(fit, pars = c("R_mu",
                                      "post_pred_c1_t1", "post_pred_c1_t2",
                                      "post_pred_c2_t1", "post_pred_c2_t2"))
# Only save random 3 subjects for plotting (to avoid RAM problems)
post_samples$post_pred_c1_t1 <- post_samples$post_pred_c1_t1[,subjs,]
post_samples$post_pred_c1_t2 <- post_samples$post_pred_c1_t2[,subjs,]
post_samples$post_pred_c2_t1 <- post_samples$post_pred_c2_t1[,subjs,]
post_samples$post_pred_c2_t2 <- post_samples$post_pred_c2_t2[,subjs,]
# Test-retest plots
# First, calculate sample test-retest Pearson's correlation (i.e. two-stage approach)
samp_data <- stan_data
samp_t1 <- rowMeans(samp_data$RT[,1,1,],na.rm=T) - rowMeans(samp_data$RT[,2,1,],na.rm=T)
samp_t2 <- rowMeans(samp_data$RT[,1,2,],na.rm=T) - rowMeans(samp_data$RT[,2,2,],na.rm=T)
samp_cor <- cor.test(samp_t1, samp_t2)
# Then Plot posterior correlation and sample test-retest
mu_delta_samples <- post_samples$R_mu[,1,2]
ggplot() +
  geom_density(aes(x = mu_delta_samples, fill = "#edafaf"), alpha = .7) +
  # scale_fill_manual(values = c("#edafaf", "#b5000c", "#700000"),
  #                  labels = c(" Normal ", "Lognormal ", "Shifted Lognormal ")) +
  geom_vline(xintercept = samp_cor$estimate, linetype = 2, linewidth = 1) +
  geom_segment(aes(x = samp_cor$conf.int[1], xend = samp_cor$conf.int[2], y = 0, yend = 0),
               color = I("black"), size = 1.5) +
  xlab("Correlation Coefficient") + ylab("density")+
  xlim(0,1) +
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "NULL",
        legend.title = element_blank())

# Posterior predictions for RTs
subj <- subjs[2]
color_scheme_set("red")
n_samples <- dim(post_samples$post_pred_c1_t1)[1]
samp <- sample(1:n_samples, 100, F)

pars1 <- post_samples$post_pred_c1_t1[samp,which(subj==subjs),]
pars2 <- post_samples$post_pred_c2_t1[samp,which(subj==subjs),]
raw1 <- samp_data$RT[subj,1,1,]
raw2 <- samp_data$RT[subj,2,1,]

p1 <- ppc_dens_overlay(y = raw1,
                       yrep = pars1) +
  geom_vline(xintercept = mean(raw1), linetype = 2, color = I("black")) +
  xlab("Congruent") + xlim(0,0.8)+ ylim(0,12)+
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none")

p2 <- ppc_dens_overlay(y = raw2,
                       yrep = pars2) +
  geom_vline(xintercept = mean(raw2), linetype = 2, color = I("black")) +
  xlab("Incongruent") + xlim(0,0.8)+ ylim(0,12)+
  theme_minimal(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.position = "none")
# p1 + p2
grid.arrange(p1, p2, ncol=2)


#########################################################################################


check_divergence <- function(fit) {

  total_iterations <- prod(dim(extract(fit, permuted=FALSE)))[1:2]
  proportion_divergent <- num_divergent / total_iterations

  message(sprintf("Number of divergent transitions: %d", num_divergent))
  message(sprintf("Proportion of divergent transitions: %f", proportion_divergent))

  if (proportion_divergent > 0.05) {
    warning("There are a significant number of divergent transitions. Consider revising your model.")
  }
}

check_divergence(fit)

# # bayesplot diagnosis
# lp_cp <- log_posterior(fit)
# head(lp_cp)
# plot(lp_cp$Value)
#
# np_cp <- nuts_params(fit)
# head(np_cp)
#
# # Extract posterior draws for later use
# posterior_cp <- as.array(fit)
#
#
# mcmc_trace(posterior_cp, pars = "R_mu", np = np_cp, window = c(300,500)) +
#   xlab("Post-warmup iteration")
#
#
# mcmc_nuts_divergence(np_cp, lp_cp)
#
# mcmc_parcoord(posterior_cp,np=np_cp)
#
# mcmc_pairs(posterior_cp, np = np_cp, pars = c("mu_mean_base","sigma_mean_base","mu_i_base[1]"),
#            off_diag_args = list(size = 0.75))

# assign to an object so we can reuse later
# scatter_theta_cp <- mcmc_scatter(
#   posterior_cp,
#   pars = c("mu_i_base[1]", "mu_sd_base"),
#   transform = list(tau = "log"), # can abbrev. 'transformations'
#   np = np_cp,
#   size = 1
# )
# scatter_theta_cp
#' Parameter Recovery for Model Testing without Experimental Data
#'
#' This function extracts DDM parameters from a Stanfit object run on simulated data. Then,
#' it and compares those parameters fitted on simulated data with the parameters used to simulate
#' the data in the first place. This function is intended to check if the model is capable to
#' "recover" parameters which are known. Plots for visual inspection are made automatically.
#' Automatic fit diagnostics yet to be implementd with this function. The comparison is done
#' for each parameter over all iterations to detect convergence problems or to spot difficult parameter
#' spaces.
#'
#' @export
#' @param fit Stanfit object fitted on synthetic data.
#' @param n_subjects Number of "Simulated Subjects".
#' @param auth_params "Ground-truth" parameters used to simulate the data. In the form of the
#' output from the function \code{\link{makeFakeParams}}.
#' @param model_name Optional name for the model to distinguish the plots and created
#' folder easier. Default is 'NAME_UNDEFINED'.
#' @return Creates a folder and saves plots which compare recovered parameters and "ground-truth" parameters
#' for each iteration of the MCMC sampler.
param_recovery <- function(fit, n_subjects, auth_params, model_name='NAME_UNDEFINED'){
  # plot divergencies by parameter:
  mat <- rstan::extract(fit, permuted=FALSE)

  n_iter <- dim(mat)[1]
  n_chains <- dim(mat)[2]

  param <- rstan::extract(fit, permuted=TRUE, inc_warmup=FALSE) %>% map(unname)

  param <- param[as.vector(map(param, length)==(n_iter*n_subjects*n_chains))]
  #filter elements in list with the amount of data only subject-level params have

  param %<>% map(function(x){as.data.frame(x)}) %>%
    map(function(x){ apply(x, 2, function(y) cumsum(y)/seq_along(y))}) %>%
    map(unname) %>%
    map(function(x){data.table::melt(x)}) %>%
    map(function(x){x[,c('Var2', 'value')]}) %>%
    map(function(x){x$Var2 <- as.factor(x$Var2); return(x)}) %>%
    melt()

  param %<>% select('Var2', 'value', 'L1')

  param$L1 <- as.factor(param$L1)

  names(param) <- c('suj', 'value', 'param')

  ##auth param
  tmpa <- c('a', 'v', 't', 'z', 'st', 'sz', 'sv')
  tmpb <- c('alpha_mu', 'delta_mu', 'tau_mu', 'beta_mu', 'tau_sigma', 'beta_sigma', 'delta_sigma')
  for(.. in 1:length(tmpa)){
    if(tmpa[..] %in% names(auth_params)){
      names(auth_params)[names(auth_params) == tmpa[..]] <- tmpb[..]
    }
  }

  auth_params$suj <- as.factor(seq.int(nrow(auth_params)))

  auth_params %<>% data.table::melt(id.vars='suj')
  names(auth_params) <- c('suj', 'param', 'value')
  param %<>% left_join(., auth_params , by=c('suj', 'param') )

  dd <- param %>% filter(param!='log_lik' & param!="lp__") %>% group_by(param) %>%  nest() %>%
    mutate(plot=map2(data, as.list(.$param), ~ggplot(data = .x ) +
                       geom_line(aes(y=value.x, x=1:length(.x$value.x))) +
                       geom_line(aes(y=value.y, x=1:length(.x$value.x)), linetype = "dashed", color='red') +
                       theme(legend.position="bottom") +
                       scale_colour_hue(name="Data Types:")  +
                       theme(legend.position="bottom") + facet_wrap(~suj, scales='free') +
                       labs(title = paste('Parameter:', .y, sep=' '),
                            subtitle = "Parameter recovery per simulated subject",
                            x = "Iterations",
                            y = 'Parameter Values',
                            caption = 'Note: Warmup iterations not included. The dashed line represents the
                                  parmeter value which generated the data.')))

  walk2(.x = dd$plot, .y = dd$param,
        ~ ggsave(device = 'pdf',
                 scale = 1, width = 12, height = 8, units = c("in"),
                 filename = paste0(model_name, '/', model_name, '_param_recovery_', .y, ".pdf"),
                 plot = .x))

}


#' Model Comparison via LOO-CV, WAIC and Raw LPPD
#'
#' Applies log-likelihood based model comparison to any number of stanfit objects and extracts LOO-CV,
#' WAIC and Raw LPPD measures and plots the models on a absolute scale ranging from 0 (best models) to
#' 1 (worst models). WAIC and LOO-CV code has been adapted directly from the "LOO" package. Because
#' most methods to prevent log-exp calculation-underflow failed, this function imputes underflow values
#' with the mean of sampled values. When this happens, the function will report the amount of imputed values
#' to inform the user. More than 5% of imputed values should not be accepted.
#' Mind how Stan saves log-likelihood values! Consult the Stan manual to check how they are saved correctly!
#' Because all models need to be loaded into memory, be wary if loading to big stanfits. Lighten the Stanfit
#' objects by discarding unnecesary iterations is advised.
#'
#' @export
#' @param stanfits At least 2 stanfit objects.
#' @param lik_name Name under which the log likelihoods have been saved in the models. Needs to be identical
#' across all Stanfit objects.
#' @param impute_inf A boolean which regulates if underflow values should be automatically imputed or not. If
#' \code{FALSE}, the models with such values will just be ignored. If \code{TRUE}, a report will be generated
#' on how many values were imputed and for which models.
#' @examples
#'
#' load(fit1.Rdata)
#' load(fit2.Rdata)
#'
#' model_judgement(fit1, fit2, impute_inf = TRUE)
#' @return Prints a table and generates a plot with the model ranked according to LOO-CV, WAIC and Raw LPPD.
model_judgement <- function(..., lik_name = "log_lik", impute_inf = TRUE) {
  #compare all models in environment: compare <- ls(pattern = 'StanDDM')
  if (!require("data.table")) {install.packages("data.table", dependencies = TRUE); library(data.table)}
  if (!require("rstan")) {install.packages("rstan", dependencies = TRUE); library(rstan)}

  arg <- list(...)
  nms <- substitute(list(...))[-1]
  names(arg) <- as.list(nms)

  if(impute_inf){
    imps <- arg}
  else{
    pre <- precheck(arg, lik_name)
    if(pre[[1]]==1){
      cat('\nAll models are OK to go.\n')
    }
    else{
      cat('\nWarning: The following models feature Log Exp underflow NaNs:\n\n')
      print(pre$troublemakers)
      cat('___________________________________\n\n')
    }
  }




  for(.. in 1:length(arg)){

    current_model <- nms[..]

    judgement <- (waic(arg[[..]], current_model, lik_name, impute_inf))

    ic <- judgement['ic']
    arg[[..]] <- list("lppd" = as.numeric(ic$ic$total['lpd']),
                      "PSIS-LOO" = as.numeric(ic$ic$total['elpd_loo']),
                      "WAIC" = as.numeric(ic$ic$total['waic']))

    if(impute_inf){
      imp <- judgement['imps']
      imps[[..]] <- list('Model Name' = as.character(current_model),
                         "%Likelihood imp." = imp$imps$lik_imp,
                         "% LOO CV imp." = imp$imps$loo_imp)
    }


  }
  if(impute_inf){
    cat('\n\n___________________________________')
    cat('\n Warning: The following models had -Inf values imputed (shows % of data imputed):\n\n')
    print(rbindlist(imps, fill=TRUE))
    cat('___________________________________\n\n')
  }

  #plotting:
  gc()
  df <- rbindlist(arg, fill=TRUE)

  madmax <- apply(abs(df), 2, FUN = max, na.rm=TRUE)

  df <- sweep(abs(df), 2, madmax, "/")
  df <- cbind(Model = names(arg), df)
  molten_df <- melt(df)
  # pd <- position_dodge(0.05)  #The errorbars overlapped, so use position_dodge to move them horizontally
  plot <- ggplot(molten_df, aes(x=variable, y=value, colour=Model,group=Model)) +
    # geom_errorbar(aes(ymin=value-se, ymax=mean+se), width=.1,position=pd) +
    geom_line() +
    geom_point()+
    geom_text(aes(label=Model),hjust=0, vjust=0,
              position = position_dodge(width=0.3),  size=2.5) +
    ylab("Absolute ranking of model fit") +
    xlab("IC")+
    scale_colour_hue(name="Models:") +
    # ylim(0, 1) +
    theme(legend.position="bottom")

  gc()
  return(list(df, plot))


}

#' LOO-CV, WAIC and Raw LPPD Calculations
#'
#' Implements LOO-CV, WAIC and Raw LPPD for the \code{\link{model_judgement}} function. Contains
#' the helper function "colVars" which calculates row-wise variances efficiently.
#'
#' @param stanfit A Stanfit object fitted on synthetic data.
#' @param current_model Name of the current model (Currently not in use.)
#' @param lik_name Name under which the log likelihoods have been saved in the models. Needs to be identical
#' across all Stanfit objects.
#' @param impute_inf A boolean which regulates if underflow values should be automatically imputed or not. If
#' \code{FALSE}, the models with such values will just be ignored. If \code{TRUE}, a report will be generated
#' on how many values were imputed and for which models.
#' @return A list with LOO-CV, WAIC and Raw LPPD calculations.
waic <- function(stanfit, current_model, lik_name, impute_inf){
  #http://kylehardman.com/BlogPosts/View/6 DIC code also from Gelman
  #Modified code of www.stat.columbia.edu/~gelman/research/unpublished/waic_stan.pdf
  #from gist.github.com/ihrke for underflow probs

  colVars <- function(a) {
    n <- dim(a)[[1]];
    c <- dim(a)[[2]];
    result <- (.colMeans(((a - matrix(.colMeans(a, n, c),
                                      nrow = n, ncol = c, byrow = TRUE)) ^ 2), n, c) * n / (n - 1))
    return(result)
  }

  log_lik <- rstan::extract(stanfit, lik_name)$log_lik


  if(impute_inf){
    lik_imp <- sum(is.infinite(log_lik))/length(log_lik)
    log_lik[is.infinite(log_lik)] <- mean(log_lik[!is.infinite(log_lik)])
  }
  dim(log_lik) <- if (length(dim(log_lik))==1) c(length(log_lik),1) else
    c(dim(log_lik)[1], prod(dim(log_lik)[2:length(dim(log_lik))]))
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  #log pointwise
  lpd <- log(colMeans(exp(log_lik))) #only when posterior simulations in Stan are correctly made (and possible)

  #waic
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic

  #loo
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))

  if(impute_inf){
    loo_imp <- sum(is.infinite(loo_weights_raw))/length(loo_weights_raw)
    loo_weights_raw[is.infinite(loo_weights_raw)] <- mean(loo_weights_raw[!is.infinite(loo_weights_raw)])
  }

  loo_weights_normalized <- loo_weights_raw/matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin (loo_weights_normalized, sqrt(S))

  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized)/colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  elpd_loo <- elpd_loo*-2
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  ic <- list(waic=total["waic"], elpd_waic=total["elpd_waic"],
             p_waic=total["p_waic"], elpd_loo=total["elpd_loo"], p_loo=total["p_loo"],
             pointwise=pointwise, total=total, se=se)

  if(impute_inf){
    imps <- list('lik_imp' = lik_imp, 'loo_imp' = loo_imp)
    return(list('ic' = ic, 'imps' = imps))
  } else{
    return(list('ic' = ic))
  }
}
