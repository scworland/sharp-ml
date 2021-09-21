library(lightgbm)
library(tidyverse)
library(tidymodels)
library(finetune)
library(treesnip)
library(rules)
library(glue)
library(jsonlite)
library(forcats)
library(butcher)
library(lobstr)
library(doMC)
library(furrr)
tidymodels_prefer()
conflicted::conflict_prefer("flatten", "purrr")
conflicted::conflict_prefer("max_rules", "rules")


#' sample diameters from log normal distribution
#' @param n integer: number of samples
#' @param beta double: meanlog = log(beta) 
#' @param alpha double: sdlog = alpha
sample_diameters <- function(n=1e6,beta=115.4,alpha=0.36){
  
  samples <- rlnorm(n, meanlog = log(beta), sdlog = alpha)
  
  return(samples)
}

#' create buckets for diameters
#' @param diameter double: groups diameters into buckets
diameter_bucket <- function(diameter){
  
  # custom case_when function for factors
  fct_case_when <- function(...) {
    args <- as.list(match.call())
    levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
    levels <- levels[!is.na(levels)]
    factor(dplyr::case_when(...), levels=levels)
  }
  
  fct_case_when(
    diameter >= 0 & diameter <= 50 ~ '<50',
    diameter > 50 & diameter <= 100 ~ '51-100',
    diameter > 100 & diameter <= 150 ~ '101-150',
    diameter > 150 & diameter <= 200 ~ '151-200',
    diameter > 200 & diameter <= 250 ~ '201-250',
    diameter > 250 & diameter <= 300 ~ '251-300',
    diameter > 300 & diameter <= 350 ~ '301-350',
    diameter > 350 ~ '>350'
  )
  
}

#' Loads csv data and does basic preprocessing steps
#' @param none
load_data <- function(){
  
  annulus <- read_csv('input_data/AnnulusData_7Jul21.csv')
  cylinder <- read_csv('input_data/CylinderData.csv')
  slab <- read_csv('input_data/Slab1Data.csv')
  cylinder_new <- read_csv('input_data/CylinderData_3Sep21.csv')
  slab_new <- read_csv('input_data/Slab1Data_3Sep21.csv')
  all_geoms_new <- read_csv('input_data/IsletData_9Sep21.csv')
  
  model_data <- bind_rows(annulus,cylinder,cylinder_new,slab,slab_new,all_geoms_new) %>%
    rename_all(tolower) %>%
    mutate(diameter = diameter*1e6,
           bucket = diameter_bucket(diameter),
           tau = ifelse(geometry == 'Annulus', tau*2, tau))
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

#' summarize k values
#' @param df: raw data
format_data <- function(raw_data) {
  
  # average diameter for each bucket
  bucket_avg <- raw_data %>%
    group_by(bucket) %>%
    summarize(diameter = mean(diameter))
  
  # average s and k for each geometry, rho, 
  # tau, and bucket. + average diameter for bucket
  formatted_data <- raw_data %>%
    group_by(geometry,bucket,tau,rho) %>%
    summarize(s = mean(s),
              k = mean(k),
              .groups='drop') %>%
    left_join(bucket_avg,by='bucket')
}

#' build workflow set of models
#' @param model_rec recipe: preprocessing steps
build_workflow_set <- function(model_rec){
  
  # linear model spec
  lm_spec <- 
    linear_reg() %>%
    set_engine('lm') %>%
    set_mode('regression')
  
  # lightGBM spec
  gbm_spec <-
    boost_tree(mtry = tune(), 
               trees = tune(), 
               min_n = tune(), 
               tree_depth = tune(),
               loss_reduction = tune(), 
               learn_rate = tune(), 
               sample_size = 0.75) %>%
    set_engine("lightgbm") %>%
    set_mode("regression")
  
  
  # xgboost spec
  xgb_spec <- 
    boost_tree(mtry = tune(), 
               trees = tune(), 
               min_n = tune(), 
               tree_depth = tune(),
               loss_reduction = tune(), 
               learn_rate = tune(), 
               sample_size = 0.75) %>%
    set_engine("xgboost") %>%
    set_mode("regression")
  
  # cubist spec
  cubist_spec <-
    cubist_rules(committees = tune(), 
                 neighbors = tune(),
                 max_rules = tune()) %>%
    set_engine('Cubist') %>%
    set_mode("regression")
  
  # neural net spec
  nnet_spec <-
    mlp(
      hidden_units = tune(),
      penalty = tune(),
      epochs = 300,
      activation = 'relu') %>%
    set_engine('nnet') %>%
    set_mode("regression")
  
  # knn spec
  knn_spec <- 
    nearest_neighbor(
      neighbors = tune(),
      weight_func = tune(),
      dist_power = tune()) %>% 
    set_engine('kknn') %>%
    set_mode("regression")
  
  ## update default parameter ranges ----
  
  # update cubist params
  cubist_params <- 
    cubist_spec %>%
    parameters() %>%
    update(
      neighbors = neighbors(range=c(0,5)),
      max_rules = max_rules(c(1,500))
    )
  
  # update gbm params
  gbm_params <- 
    gbm_spec %>%
    parameters() %>%
    update(
      mtry = mtry(range(c(1,11))),
      trees = trees(range=c(50,1500)),
      min_n = min_n(range=c(5,100))
    )
  
  # update xgb params
  xgb_params <- 
    xgb_spec %>%
    parameters() %>%
    update(
      mtry = mtry(range(c(1,11))),
      trees = trees(range=c(50,1500)),
      min_n = min_n(range=c(5,100))
    )
  
  # update nnet params
  nnet_params <- 
    nnet_spec %>%
    parameters() %>%
    update(
      hidden_units = hidden_units(range=c(5,25))
    )
  
  # update knn params
  knn_params <- 
    knn_spec %>%
    parameters() %>%
    update(
      weight_func = weight_func(values=c('gaussian','rectangular'))
    )
  
  ## Create workflow set ----
  
  # model spec list
  model_specs <- list(
    lm = lm_spec,
    nnet = nnet_spec,
    knn = knn_spec,
    xgb = xgb_spec,
    gbm = gbm_spec,
    cubist = cubist_spec
  )
  
  # workflow set 
  wkflow_set <- workflow_set(
    preproc = list(base_rec = model_rec), 
    models = model_specs,
    cross = TRUE) %>%
    mutate(wflow_id = names(model_specs)) %>%
    option_add(param_info = cubist_params, id = "cubist") %>%
    option_add(param_info = gbm_params, id = "gbm") %>%
    option_add(param_info = nnet_params, id = 'nnet') %>%
    option_add(param_info = xgb_params, id = 'xgb') %>%
    option_add(param_info = knn_params, id = 'knn')
  
  return(wkflow_set)
  
}

#' tune the models using grid sampling
#' @param workflow_set workflow_set: workflow set of models
#' @param cv_splits vfold_cv: cross validated folds
#' @param n integer: number of hyperparemerters for grid sampling
tune_models <- function(workflow_set,cv_splits,n=50){
  
  # tuning controls
  tune_ctrl <-
    control_grid(
      save_pred = FALSE,
      save_workflow = FALSE,
      allow_par = TRUE,
      parallel_over = 'everything'
    )
  
  # train models and store prediction metrics
  options(tidymodels.dark = TRUE)
  tune_results <- workflow_set %>%
    workflow_map(
      fn="tune_grid",
      seed = 123,
      resamples = cv_splits,
      grid = n,
      control = tune_ctrl,
      verbose=TRUE
    )
  
  return(tune_results)
}

#' workaround for: https://github.com/curso-r/treesnip/issues/33
#' @param model workflow: final lightgbm workflow
#' @geometry string: geometry of device
save_lightgbm <- function(model,geometry){
  
  # save workflow
  saveRDS(model, glue('r_objects/production/{geometry}_gbm_fit_workflow'))
  
  # extract lightgbm model object
  lightgbm_fit = extract_fit_parsnip(model)
  
  # save model object
  lgb.save(lightgbm_fit$fit, glue('r_objects/production/{geometry}_gbm_fit_model'))
}

#' fit on resamples using best combination of params
#' @param model_id string: string to identify model
#' @param tune_results workflow_set: results of tuning grid
#' @param cv_splits vfold_cv: cross validation folds
fit_workflow_on_resamples <- function(model_id,tune_results,cv_splits,formatted_data){
  
  # extract meta data
  meta_data <- cv_splits$splits[1][[1]]$data %>%
    select(-s)
  
  # geometry
  geometry <- unique(meta_data$geometry)
  
  # final hyper parameters
  model_params <- tune_results %>% 
    extract_workflow_set_result(model_id) %>% 
    select_best(metric = "rmse")
  
  # final model
  model_final <- 
    tune_results %>% 
    extract_workflow(model_id) %>% 
    finalize_workflow(model_params) 
  
  # fit model to formatted data and save
  final_fit <- model_final %>% 
    fit(formatted_data)
  
  if(model_id == 'gbm'){
    final_fit %>% 
      save_lightgbm(geometry)
  }else{
    final_fit %>% 
      saveRDS(glue('r_objects/production/{geometry}_{model_id}_fit'))
  }
  
  # save predictions
  cntrl <- control_resamples(verbose=TRUE,save_pred=TRUE)
  
  # resample fits
  resample_fits <- model_final %>%
    fit_resamples(cv_splits,control=cntrl)
  
  # cv errors
  cv_error <- resample_fits %>%
    collect_metrics() %>%
    mutate(model_id=model_id,
           geometry=geometry) %>%
    select(geometry, model_id,metric=.metric,mean,std_err,folds=n)
  
  # extract predictions
  cv_predictions <- resample_fits %>%
    collect_predictions(summarize=TRUE) %>%
    mutate(model_id = model_id) %>%
    bind_cols(meta_data) %>%
    mutate(s_est = ifelse(.pred>1,1,(ifelse(.pred<0,0,.pred))),
           k_est = s_est*(diameter*1e-6)^3/((150e-6)^3)) %>%
    select(model_id, geometry, bucket, rho, tau, diameter,
           s_est, s, k_est, k)
  
  # model to export parameters
  model_params <- model_params %>%
    mutate(model_id = model_id,
           geometry=geometry)
  
  results <- list(cv_error = cv_error,
                  cv_predictions = cv_predictions,
                  model_params = model_params)
  
  return(results)
  
}

#' run model pipeline
#' @param geometry_subset string: geometry for pipeline
#' @param raw_data df: raw data
#' @param model_rec recipe: preprocessing template
#' @param v integer: number of folds for cross validation
#' @param grid_n integer: number of points to sample for grid
model_train_pipeline <- function(geometry_subset, raw_data, model_rec, v=10, grid_n=50){
  
  print(glue::glue('Training models for {geometry_subset} geometry'))
  
  # formatted data
  formatted_data <- format_data(raw_data) %>%
    filter(geometry == geometry_subset)
  
  # cross validation
  cv_splits <- vfold_cv(formatted_data, v = v)
  
  # build workflow set
  workflow_set <- build_workflow_set(model_rec)
  
  # tune the models using CV
  tune_results <- tune_models(workflow_set,cv_splits,n=grid_n)
  
  # gather cv predictions using best params for each model
  resample_fits <- workflow_set$wflow_id %>% 
    map(fit_workflow_on_resamples,tune_results,cv_splits,formatted_data)
  
}

#' calculate model weights
#' @param cv_predictions df: dataframe of cross validated predictions
#' @save boolean: save model weights
calculate_model_weights <- function(cv_predictions,save=TRUE){
  
  model_weights <- cv_predictions %>%
    group_by(geometry) %>%
    mutate(tau_bucket = ifelse(tau <= median(tau),glue('tau <= {median(tau)}'),glue('tau > {median(tau)}')),
           rho_bucket = ifelse(rho <= median(rho),glue('rho <= {median(rho)}'),glue('rho > {median(rho)}'))) %>%
    ungroup() %>%
    group_by(model_id,geometry,bucket,rho_bucket,tau_bucket) %>%
    summarize(r2 = abs(cor(s,s_est)),
              .groups='drop') 
  
  if(save){
    saveRDS(model_weights,'r_objects/production/model_weights')
  }
  
  return(model_weights)
}

#' plot CV predictions
#' @param cv_predictions df: dataframe of cross validated predictions
plot_predictions <- function(cv_predictions) {
  
  cv_predictions %>%
    group_by(geometry,bucket,model_id) %>%
    summarize(med_actual = mean(k),
              med_estimate = mean(k_est),
              .groups='drop') %>%
    pivot_longer(cols = contains(c('med_')),
                 values_to = 'med',
                 names_to='type',
                 names_prefix='med_') %>%
    left_join(
      cv_predictions %>%
        group_by(geometry,bucket,model_id) %>%
        summarize(sd_actual = sd(k),
                  sd_estimate = sd(k_est),
                  .groups='drop') %>%
        pivot_longer(cols = contains(c('sd_')),
                     values_to='sd',
                     names_to='type',
                     names_prefix='sd_'),
      by=c('geometry','bucket','type','model_id')
    ) %>%
    mutate(geometry = ifelse(geometry=='Slab1','Planar slab',geometry)) %>%
    ggplot() +
    geom_bar(aes(x=bucket,y=med,fill=type),
             stat='identity',position=position_dodge(), color='black') +
    facet_grid(model_id~geometry) +
    scale_fill_manual(values=c('#F7F7F7','#B31B1B'), labels=c('SHARP','SHARP-ML')) +
    geom_errorbar(aes(x=bucket,ymin=med-0.5*sd, ymax=med+0.5*sd, 
                      group=type), 
                  width=0.2,position=position_dodge(0.9)) +
    labs(y = expression(kappa), fill = '', 
         x=expression('Islet Size Group ('*mu*'m)')) +
    theme_minimal(base_size=14)
}


