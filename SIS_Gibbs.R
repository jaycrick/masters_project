library(pacman)
pacman::p_load("deSolve", "tidyverse", "latex2exp")

# Task 1 ----
set.seed(4503)
SIS_gamma <- 1/rgamma(1, 12, 4) # 0.3576638
SIS_beta <- rgamma(1, 3, 4/3/SIS_gamma) # 1.146848
R_0 = SIS_beta/SIS_gamma; R_0 # 3.206499

R_0_obs = rpois(1, R_0); R_0_obs

no_samples = 5000
init_beta = 1
init_gamma = 0.3
init_sample = c(beta_samp = init_beta, gamma_samp = init_gamma)

gibbs_samples = as_tibble_row(init_sample)

for (i in 2:no_samples){
  new_gamma = 1/rgamma(1, 12 + R_0_obs, 4 + gibbs_samples$beta_samp[i-1])
  new_beta = rgamma(1, 3 + R_0_obs, 4/3/SIS_gamma + 1/new_gamma)
  gibbs_samples <- add_row(gibbs_samples, beta_samp = new_beta, gamma_samp = new_gamma)
}

gibbs_samples = gibbs_samples %>% mutate(R0 = beta_samp/gamma_samp)

ggplot(gibbs_samples, mapping = aes(x = beta_samp, y = gamma_samp)) +
  geom_point(alpha = 0.1) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", colour="white") +
  theme_bw() +
  scale_x_continuous(expand = expansion(0), limits = c(0, 3)) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 1.5)) +
  geom_path(data = head(gibbs_samples, 15), colour = "red") +
  geom_point(data = head(gibbs_samples, 15), colour = "red") +
  labs(
    x = TeX(r"(\beta)"),
    y = TeX(r"(\gamma)")
  )
  

ggsave("C:/Users/jckricket/Dropbox/Apps/Overleaf/M_Scimat_Thesis/images/SIS_gibbs.pdf", width = 5, height = 5)
