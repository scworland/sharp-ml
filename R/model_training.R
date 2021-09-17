## Setup ----

# source utility functions
source('production_scripts/training_utils.r')

# register cores for parallel processing
registerDoMC(cores = 6)

# load data
raw_data <- load_data() 

# plot training data
plot_heatmap(raw_data, s)

# create template for preprocessing
template_data <- raw_data %>%
  format_data() %>%
  head()

# preprocessing template
model_rec <-
  recipe(s ~ tau + rho + diameter + bucket + geometry, data = template_data) %>%
  update_role(bucket, new_role = "bucket") %>%
  update_role(geometry, new_role = "geometry") %>%
  step_normalize(all_numeric_predictors()) 

# tune, train, and test models using cross validation
resampled_fits <- unique(raw_data$geometry) %>%
  map(model_train_pipeline, raw_data, model_rec, v=10, grid_n=100)

# cross validation metrics 
cv_metrics <- resampled_fits %>%
  map_dfr(~map(.x, ~.x$cv_error))

# cross validation predictions
cv_predictions <- resampled_fits %>% 
  map_dfr(~map(.x, ~.x$cv_predictions))

# model parameters
model_parameters <- resampled_fits %>% 
  map(~map(.x, ~.x$model_params))

# plot cross validation predictions
cv_predictions %>%
  plot_predictions()

# calculate model weights for ensemble
model_weights <- calculate_model_weights(cv_predictions,save=TRUE)


