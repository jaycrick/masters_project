library(pacman)
pacman::p_load("deSolve", "tidyverse")

# Task 1 ----
set.seed(4503)
SIS_gamma <- 1/rgamma(1, 12, 4) # 0.3596567
SIS_beta <- rgamma(1, 3, 4/3/SIS_gamma) # 1.146848
R_0 = SIS_beta/SIS_gamma; R_0 # 3.206499

R_0_obs = rpois(1, R_0); R_0_obs

no_samples = 50000
init_beta = 1
init_gamma = 0.3
init_sample = c(beta_samp = init_beta, gamma_samp = init_gamma)

abc_samples = tribble(~beta_samp, ~gamma_samp)

while (dim(abc_samples)[1] < no_samples){
  new_gamma = 1/rgamma(1, 12, 4)
  new_beta = rgamma(1, 3, 4/3/SIS_gamma)
  R0_guess = rpois(1, new_beta/new_gamma)
  if (R0_guess == R_0_obs){
    abc_samples <- add_row(abc_samples, beta_samp = new_beta, gamma_samp = new_gamma)
  }
}

abc_samples = abc_samples %>% mutate(R0 = beta_samp/gamma_samp)
mean(abc_samples$R0)

ggplot(abc_samples, mapping = aes(x = beta_samp, y = gamma_samp)) +
  geom_point(alpha = 0.1) +
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  theme_bw() +
  scale_x_continuous(expand = expansion(0), limits = c(0, 3)) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 1.5))

# ggplot(abc_samples, mapping = aes(x = gamma_samp)) +
#   geom_histogram(bins = 80) +
#   theme_bw() +
#   xlim(c(0, 1))

ggsave("C:/Users/jckricket/Dropbox/Apps/Overleaf/M_Scimat_Thesis/images/SIS_ABC.pdf", width = 5, height = 5)
