# load results-----------
load('./sardine_study/results/app_evidence_power_posterior_M2.RData')

app_evidence_power_posterior_M2 <- app_evidence_power_posterior_M3
rm(app_evidence_power_posterior_M3)

load('./sardine_study/results/app_evidence_power_posterior_M3.RData')

app_evidence_power_posterior_M3 <- app_evidence_power_posterior_M1
rm(app_evidence_power_posterior_M1)




load('./sardine_study/results/app_evidence_power_posterior_M1.RData')
app_evidence_power_posterior_M1 <-
  app_evidence_power_posterior_mixed_linear_model_global_sardine_fishing
rm(app_evidence_power_posterior_mixed_linear_model_global_sardine_fishing)

# summary------------

mean(app_evidence_power_posterior_M1); sd(app_evidence_power_posterior_M1)
mean(app_evidence_power_posterior_M2); sd(app_evidence_power_posterior_M2)
mean(app_evidence_power_posterior_M3); sd(app_evidence_power_posterior_M3)



