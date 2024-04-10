library("deSolve")
library("tidyverse")
library("viridis")

pv_champ_function <- function(t, state, parameters) { # with optional birth death
  with(as.list(c(state, parameters)), {
    # rate of change
    ds_0 <- -(1 - alpha * beta) * (lambda * (i_L + i_0) + delta) * s_0 + 
      (lambda * (i_0 + i_L) + delta) * alpha * beta * s_L + 
      alpha*beta*f*s_L + 
      gamma_L*s_L + 
      r*i_0
    
    di_0 <- -(lambda * (i_0 + i_L) + delta) * i_0 + 
      gamma_L*i_L - r*i_0
    
    ds_L <- -(1 - alpha*(1 - beta))*(lambda*(i_L + i_0) + delta + f)*s_L +
      alpha*(1-beta)*(lambda*(i_0+i_L) + delta)*s_0 - gamma_L*s_L + r*i_L
    
    di_L <- (1 - alpha)*(lambda*(i_L + i_0) + delta)*(s_0 + s_L) +
      (lambda*(i_L + i_0) + delta)*i_0 + (1 - alpha)*f*s_L - gamma_L*i_L - r*i_L
    
    # return the rate of change
    list(c(ds_0, di_0, ds_L, di_L))
  }) # end with(as.list ...
}

pv_champ_alpha <- 0.4 #prop of effective care
pv_champ_beta <- 0.4 #prop of radical cure
pv_champ_gamma_L <- 1 / 223 #liver stage clearance rate
pv_champ_delta <- 0.05 #prop of imported cases
pv_champ_lambda <- 0.04 #transmission rate
pv_champ_f <- 1 / 72 #relapse frequency
pv_champ_r <- 1 / 60 #blood stage clearance rate

pv_champ_parameters <- c(alpha = pv_champ_alpha, beta = pv_champ_beta, gamma_L = pv_champ_gamma_L, delta = pv_champ_delta, lambda = pv_champ_lambda, f = pv_champ_f, r = pv_champ_r)
state <- c(s_0 = 0.99, i_0 = 0.01, s_L = 0, i_L = 0)

pand_length <- 40

times <- seq(0, pand_length, by = 2)

pv_champ_outputs <- ode(y = state, times = times, func = pv_champ_function, parms = pv_champ_parameters) %>%
  as_tibble() %>%
  mutate(across(.cols = everything(), .fns = as.numeric))


pv_champ_long_outputs <- pivot_longer(pv_champ_outputs, -time, names_to = "state", values_to = "proportions")

pv_champ_plot <- ggplot(pv_champ_long_outputs, aes(x = time, y = proportions, colour = state)) +
  geom_line(size = 0.75) +
  theme_bw() +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since first infection") +
  ylab("Proportion in each compartment")
pv_champ_plot