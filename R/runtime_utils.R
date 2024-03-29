library(tidyverse)
library(tidymodels)
library(treesnip)
library(rules)
library(glue)
tidymodels_prefer()
conflicted::conflict_prefer("flatten", "purrr")
conflicted::conflict_prefer("max_rules", "rules")

#' create prediction set
#' @param geometry string: geometry of device
#' @param rho double: rho for prediction
#' @param tau double: tau for prediction
create_prediction_set <- function(geometry,rho,tau){
  
  # load average diameters from 1e6 sample
  bucket_average <- readRDS('r_objects/production/bucket_average')
  
  prediction_data <- tibble(
    geometry = geometry,
    rho = rho/100,
    tau = tau/1000
  ) %>%
    expand_grid(bucket_average) 
  
}

#' workaround for: https://github.com/curso-r/treesnip/issues/33
#' @geometry string: geometry of device
load_lightgbm <- function(geometry){
  
  model <- readRDS(glue('r_objects/production/{geometry}_gbm_fit_workflow'))
  model_lgb <- lightgbm::lgb.load(glue('r_objects/production/{geometry}_gbm_fit_model'))
  model$fit$fit$fit <- model_lgb
  
  return(model)
}

#' load the models
#' @param geometry string: geometry of device
load_models <- function(geometry){
  
  # load geometry specific models
  cubist <- readRDS(glue('r_objects/production/{geometry}_cubist_fit'))
  knn <- readRDS(glue('r_objects/production/{geometry}_knn_fit'))
  lm <- readRDS(glue('r_objects/production/{geometry}_lm_fit'))
  nnet <- readRDS(glue('r_objects/production/{geometry}_nnet_fit'))
  xgb <- readRDS(glue('r_objects/production/{geometry}_xgb_fit'))
  gbm <- load_lightgbm(geometry)
  
  model_list <- list(cubist=cubist,knn=knn,lm=lm,nnet=nnet,gbm=gbm,xgb=xgb)
  
  return(model_list)
}

#' make predictions
#' @param geometry string: geometry of device
#' @param rho double: rho for prediction
#' @param tau double: tau for prediction
#' @param ensemble boolean: should predictions be ensembled
make_predictions <- function(geometry,rho,tau,ensemble=TRUE){
  
  # create prediction data
  prediction_data <- create_prediction_set(geometry,rho,tau)
  
  # load models
  models <- load_models(geometry)
  
  # load model weights
  model_weights <- readRDS('r_objects/production/model_weights')
  
  geometry_medians <- model_weights %>%
    mutate(rho_median = parse_number(rho_bucket),
           tau_median = parse_number(tau_bucket)) %>%
    distinct(geometry, rho_median, tau_median)
  
  
  model_prediction <- list()
  for(i in 1:length(models)){
    model_prediction[[i]] <- prediction_data %>%
      mutate(s_est = predict(models[[i]],prediction_data)$.pred,
             k_est = s_est*(diameter*1e-6)^3/((150e-6)^3),
             model_id = names(models[i])) %>%
      select(geometry,rho,tau,bucket,model_id,diameter,s_est,k_est)
  }
  
  if(ensemble){
    model_predictions <- bind_rows(model_prediction) %>%
      left_join(geometry_medians, by='geometry') %>%
      mutate(tau_bucket = ifelse(tau <= tau_median,glue('tau <= {tau_median}'),glue('tau > {tau_median}')),
             rho_bucket = ifelse(rho <= rho_median,glue('rho <= {rho_median}'),glue('rho > {rho_median}'))) %>%
      left_join(model_weights, 
                by = c('geometry','bucket','model_id','tau_bucket','rho_bucket')) %>%
      group_by(geometry,bucket,tau,rho) %>%
      mutate(w = r2/sum(r2)) %>%
      summarize(s_mu = round(sum(s_est*w),4),
                s_sigma = round(sum(w*(s_est - s_mu)^2),4),
                k_mu = round(sum(k_est*w),4),
                k_sigma = round(sum(w*(k_est - k_mu)^2),4),
                .groups='drop') %>%
      mutate(s_mu = ifelse(s_mu > 1, 1, s_mu),
             k_sigma = ifelse(k_sigma > 0.75*k_mu, 0.75*k_mu, k_sigma),
             kieq = c(NA, 0.167,0.667,1.685,3.5,6.315,10.352,15.833)) 
  }else{
    model_predictions <- bind_rows(model_prediction)
  }
  
  return(model_predictions)
  
}

#' predict over grid of rho and tau
#' @param raw_data df: needed to find range of variables
#' @param n int: number of grid points
predict_over_grid <- function(raw_data, n=15){
  
  # create rho and tau grids
  tau_grid <- raw_data %>%
    group_by(geometry) %>%
    summarize(tau_min = min(tau)*1000,
              tau_max = max(tau)*1000,
              .groups='keep') %>%
    group_modify(~ seq(from=.x$tau_min,to=.x$tau_max,length.out=n) %>%
                   tibble::enframe(value = "tau"),
                 .keep=TRUE) %>%
    select(-name) 
  
  rho_grid <- raw_data %>%
    group_by(geometry) %>%
    summarize(rho_min = min(rho)*100,
              rho_max = max(rho)*100,
              .groups='keep') %>%
    group_modify(~ seq(from=.x$rho_min,to=.x$rho_max,length.out=n) %>%
                   tibble::enframe(value = "rho"),
                 .keep=TRUE) %>%
    select(-name)
  
  expanded_grids <- tau_grid %>%
    inner_join(rho_grid,by='geometry') %>%
    mutate(ensemble = TRUE)
  
  # make predictions over the grid
  grid_predictions <- expanded_grids %>%
    pmap_dfr(make_predictions) %>% 
    group_by(geometry,rho,tau) %>%
    summarize(s_mu = mean(s_mu), 
              s_sigma = mean(s_sigma),
              k_mu = mean(k_mu),
              k_sigma = mean(k_sigma),
              .groups='drop')
  
  
  return(grid_predictions)
  
}

#' heatmap of input data
#' @param data df: input data
#' @param var df: variable to plot
plot_heatmap <- function(data,var){
  
  summary_var <- enquo(var)
  summary_nm <- as_label(summary_var)
  
  data %>%
    mutate(rho=rho,tau=tau*1000) %>%
    group_by(rho,tau,geometry) %>%
    summarize(!!summary_nm := mean({{summary_var}})) %>%
    ggplot() +
    geom_tile(aes(rho,tau,fill=!!summary_var)) +
    facet_wrap(~geometry, scales='free_y', ncol = 1) +
    scale_fill_viridis_c() +
    scale_x_continuous(labels=scales::percent) +
    theme_bw()
  
}
