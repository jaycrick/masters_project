set.seed(1094375)

library(tidyverse)
library(here)

# Coin Toss Example ----------

prob_heads = runif(1)

tosses = 10

heads = rbinom(1, tosses, prob_heads)

rho_accept = function(x, y) {
  if (y < 0 | y > 1){
    return(0)
  }
  min(dbinom(heads, tosses, y)/dbinom(heads, tosses, x), 1)
}

no_samples = 30000
init_sample = 0.5

Norm_MH_samples = c(init_sample)

for (i in 2:no_samples){
  last_sample = Norm_MH_samples[i-1]
  new_guess = rnorm(1, mean = last_sample, sd = sqrt(1/12))
  accept_prob = rho_accept(last_sample, new_guess)
  Norm_MH_samples[i] = sample(c(last_sample, new_guess), 1, prob = c(1- accept_prob, accept_prob))
}

Unif_MH_samples = c(init_sample)

for (i in 2:no_samples){
  last_sample = Unif_MH_samples[i-1]
  new_guess = runif(1, last_sample - 1/2, last_sample + 1/2)
  accept_prob = rho_accept(last_sample, new_guess)
  Unif_MH_samples[i] = sample(c(last_sample, new_guess), 1, prob = c(1- accept_prob, accept_prob))
}

MH_tbl = tibble("N(y, 1/12)" = Norm_MH_samples, "Unif(y - 1/2, y+1/2)" = Unif_MH_samples) %>%
  pivot_longer(everything(), names_to = "Proposal Distribution q(x|y)")

ggplot(MH_tbl, mapping = aes(x = value, color = `Proposal Distribution q(x|y)`)) +
  geom_function(mapping = aes(linetype = "True Posterior"), fun = dbeta, args = list(shape1 = 1 + heads, shape2 = 1 + tosses - heads), colour = "black") +
  geom_function(mapping = aes(linetype = "Prior"), fun = ~1, colour = "black", linetype = "dashed") +
  geom_freqpoly(aes(y = after_stat(density)), binwidth = 0.03, position = "dodge") +
  scale_linetype_manual(
    values = c(
      "Prior" = "dashed",
      "True Posterior" = "solid"
               )
  ) +
  theme_bw() +
  theme(text = element_text(family = "serif"), 
        legend.position = "bottom",
        legend.box="vertical",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 3)) +
  scale_x_continuous(expand = expansion(0), limits = c(0,1)) +
  labs(title = "Metropolis Hastings Coin Toss Posterior samples", x = "Probability of Heads", y = "Density", fill=NULL, linetype = NULL)

ggsave("C:/Users/jckricket/Dropbox/Apps/Overleaf/M_Scimat_Thesis/images/coin_MH_R.pdf", width = 5, height = 5)


