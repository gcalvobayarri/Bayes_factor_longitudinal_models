load("./simulation study/results/app_evidence_power_posterior_fixed_effects_model_simulated_data.RData")
load("./simulation study/results/app_evidence_power_posterior_mixed_linear_model_simulated_data.RData")
load("./simulation study/results/app_evidence_power_posterior_random_intercept_model_simulated_data.RData")
load("./simulation study/results/app_evidence_power_posterior_random_trend_model_simulated_data.RData")

df_results <- data.frame(log_evidence = c(app_evidence_power_posterior_fixed_effects_model_simulated_data,
                                          app_evidence_power_posterior_random_intercept_model_simulated_data,
                                          app_evidence_power_posterior_random_trend_model_simulated_data,
                                          app_evidence_power_posterior_mixed_linear_model_simulated_data),
                         Model = rep(1:4, each = 10))

library(ggplot2)
ggplot(df_results, aes(group=Model, y=log_evidence)) + 
  geom_boxplot()


tapply(df_results$log_evidence, df_results$Model, summary)
tapply(df_results$log_evidence, df_results$Model, mean)
tapply(df_results$log_evidence, df_results$Model, sd)
