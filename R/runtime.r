
source('production_scripts/runtime_utils.r')

# parameters
geometry = 'Cylinder'
rho = 14.5
tau = 1.3

# prediction for set of parameters
make_predictions(geometry,rho,tau,ensemble=TRUE)

# predict over grid of parameters
grid_predictions <-load_data() %>%
  predict_over_grid(n=15)

# plot heatmap of predictions
grid_predictions %>%
  plot_heatmap(s)