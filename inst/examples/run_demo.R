
# Example of BATON API
library(BATON)
sim_fun <- function(theta, fidelity=c("low","high")) {
  c(power=0.85, type1=0.08, EN=150, ET=20)
}
bounds <- list(eff=c(0.8,0.99), fut=c(0.05,0.3), ev=c(5,50), nmax=c(120,240))
res <- bo_calibrate(sim_fun, bounds, objective="EN",
                    constraints=list(power=c("ge",0.8), type1=c("le",0.1)))
print(res$best_theta)
