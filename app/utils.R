library(tidyverse)
library(tidymodels)
library(treesnip)
library(Cubist)
library(lightgbm)
library(xgboost)
library(nnet)
library(kknn)
library(rules)
library(scales)
library(glue)
library(jsonlite)
tidymodels_prefer()
conflicted::conflict_prefer("flatten", "purrr")
conflicted::conflict_prefer("max_rules", "rules")


#' define credentials
#' @param none

define_credentials <- function(){
  
  credentials <- data.frame(
    user = c("aernst", "sworland","guest"), # mandatory
    password =  c("riley robb","riley robb","sharp-ml"), # mandatory
    admin = c(FALSE, TRUE, FALSE),
    comment = "Simple and secure authentification mechanism 
  for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
  )
  
  create_db(
    credentials_data = credentials,
    sqlite_path = "credentials/database.sqlite", # will be created
    passphrase = NULL
  )
}

#' create prediction set
#' @param geometry string: geometry of device
#' @param rho double: rho for prediction
#' @param tau double: tau for prediction
create_prediction_set <- function(geometry,rho,tau){
  
  # load average diameters from 1e6 sample
  bucket_average <- readRDS('model_objects/bucket_average')
  
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
  
  model <- readRDS(glue('model_objects/{geometry}_gbm_fit_workflow'))
  model_lgb <- lightgbm::lgb.load(glue('model_objects/{geometry}_gbm_fit_model'))
  model$fit$fit$fit <- model_lgb
  
  return(model)
}

#' load the models
#' @param geometry string: geometry of device
load_models <- function(geometry){
  
  # load geometry specific models
  cubist <- readRDS(glue('model_objects/{geometry}_cubist_fit'))
  knn <- readRDS(glue('model_objects/{geometry}_knn_fit'))
  lm <- readRDS(glue('model_objects/{geometry}_lm_fit'))
  nnet <- readRDS(glue('model_objects/{geometry}_nnet_fit'))
  xgb <- readRDS(glue('model_objects/{geometry}_xgb_fit'))
  gbm <- load_lightgbm(geometry)
  
  model_list <- list(cubist=cubist,knn=knn,lm=lm,nnet=nnet,gbm=gbm,xgb=xgb)
  
  return(model_list)
}

#' make predictions
#' @param geometry string: geometry of device
#' @param rho double: rho for prediction
#' @param tau double: tau for prediction
#' @param ensemble boolean: should predictions be ensembled
make_predictions <- function(geometry,rho,tau,ensemble=TRUE,updateProgress=NULL){
  
  # create prediction data
  prediction_data <- create_prediction_set(geometry,rho,tau)
  
  # load model weights
  model_weights <- readRDS('model_objects/model_weights')
  
  geometry_medians <- model_weights %>%
    mutate(rho_median = parse_number(rho_bucket),
           tau_median = parse_number(tau_bucket)) %>%
    distinct(geometry, rho_median, tau_median)
  
  if (is.function(updateProgress)) {
    updateProgress(value = 1, detail = 'Generating predictions')
    Sys.sleep(0.5)
  }
  
  # load models
  models <- load_models(geometry)

  if (is.function(updateProgress)) {
    updateProgress(value = 2, detail = 'Generating predictions')
    Sys.sleep(1)
  }
  
  model_prediction <- list()
  for(i in 1:length(model_weights)){
    model_prediction[[i]] <- prediction_data %>%
      mutate(s_est = predict(models[[i]],prediction_data)$.pred,
             k_est = s_est*(diameter*1e-6)^3/((150e-6)^3),
             model_id = names(models[i])) %>%
      select(geometry,rho,tau,bucket,model_id,diameter,s_est,k_est)
  }
  
  if (is.function(updateProgress)) {
    updateProgress(value = 3, detail = 'Generating predictions')
    Sys.sleep(0.5)
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
  
  if (is.function(updateProgress)) {
    updateProgress(value = 4, detail = 'Modeling complete!')
  }
  
  return(model_predictions)
  
}


# plot predictions
plot_predictions <- function(formatted_predictions, geometry, rho, tau){

  geometry <- ifelse(geometry == 'Slab1','Planar slab',geometry)
  
  formatted_predictions %>%
    ggplot(aes(x=bucket)) +
    geom_errorbar(aes(ymin=k_mu - k_sigma, ymax = k_mu + k_sigma),
                  width=0.1) +
    geom_point(aes(y=k_mu), fill='white', shape = 21) +
    labs(y=expression(kappa),
         x=expression('Islet Size Group ('*mu*'m)'),
         title = glue::glue("SHARP-ML™ prediction for {geometry} geometry"),
         caption = bquote(rho~"="~.(rho)*"% and "~tau~"="~.(tau)~"mm")) +
    theme_bw(base_size=14) +
    theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
  
}


