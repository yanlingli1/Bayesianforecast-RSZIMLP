# README

This README file lists and describes the simulation materials for forecasting zero-inflated count time series with the regime-switching zero-inflated multilevel Poisson (RS-ZIMLP) model.


Reference:

Li, Y., Oravecz, Z., Zhou, S., Bodovski, Y., Barnett, I.J., Chi, G., Zhou, Y., Friedman, N.P., Vrieze, S.I., & Chow, S.-M. (2022). Bayesian Forecasting with a Regime-Switching Zero-Inflated Multilevel Poisson Regression Model: An Application to Adolescent Alcohol Use with Spatial Covariates. *Psychometrika, 87*(2), 376-402, DOI: [10.1007/s11336-021-09831-9](https://doi.org/10.1007/s11336-021-09831-9).

- *RSZIMLP_Data_Generation.R*: code for simulating data based on an RS-ZIMLP model  
	- *SimulatedData_RSZIMLP_moderate_N200T60_1.Rdata*: example simulated dataset
- *RSZIMLP_JAGS_Code.R*: JAGS code for implementing one-step forecasts on the simulated multi-subject time series data based on the RS-ZIMLP model
- *postcalc.R*: code for calculating summary statistics (e.g., means, medians, modes, standard deviations, credible intervals, effective sample sizes, and Rhat statistics) of the posterior distributions
