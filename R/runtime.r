
source('production_scripts/runtime_utils.r')

geometry = 'Cylinder'
rho = 14.5
tau = 1.3

geometry = 'Slab1'
rho = 6.5
tau = 0.65


make_predictions(geometry,rho,tau,ensemble=TRUE)

grid_predictions <- predict_over_grid(raw_data,n=15)

grid_predictions %>%
  plot_heatmap(s)

grid_predictions %>%
  group_by(rho,tau,geometry) %>%
  summarize(s = mean(s),
            .groups='drop') %>%
  ggplot(aes(rho,tau,fill=s)) +
  geom_point(shape=21) +
  #geom_tile() +
  facet_wrap(~geometry, scales='free_y', ncol = 1) +
  scale_fill_viridis_c() +
  scale_x_continuous(labels=scales::percent) +
  theme_bw()




    
    






