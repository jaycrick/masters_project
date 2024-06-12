set.seed(1094375)

library(tidyverse)
library(here)

# Coin Toss Example ----------

tosses <- 10

heads <- 6

rho_accept <- function(x, y) {
  if (y < 0 | y > 1) {
    return(0)
  }
  min(dbinom(heads, tosses, y) / dbinom(heads, tosses, x), 1)
}

no_samples <- 150000
init_sample <- 0.5

Norm_MH_samples <- c(init_sample)

for (i in 2:no_samples) {
  last_sample <- Norm_MH_samples[i - 1]
  new_guess <- rnorm(1, mean = last_sample, sd = sqrt(1 / 12))
  accept_prob <- rho_accept(last_sample, new_guess)
  Norm_MH_samples[i] <- sample(c(last_sample, new_guess), 1, prob = c(1 - accept_prob, accept_prob))
}

Unif_MH_samples <- c(init_sample)

for (i in 2:no_samples) {
  last_sample <- Unif_MH_samples[i - 1]
  new_guess <- runif(1, last_sample - 1 / 2, last_sample + 1 / 2)
  accept_prob <- rho_accept(last_sample, new_guess)
  Unif_MH_samples[i] <- sample(
    c(last_sample, new_guess), 
    1, 
    prob = c(1 - accept_prob, accept_prob))
}

MH_tbl <- tibble(
  "N(p, 1/12)" = Norm_MH_samples[seq(1001, no_samples, 5)],
  "Unif(p - 1/2, p+1/2)" = Unif_MH_samples[seq(1001, no_samples, 5)]
) %>%
  pivot_longer(everything(), names_to = "Proposal\nDistribution q")

ggplot(MH_tbl, mapping = aes(x = value, color = `Proposal\nDistribution q`)) +
  geom_function(
    mapping = aes(linetype = "True Posterior"),
    fun = dbeta,
    args = list(shape1 = 1 + heads, shape2 = 1 + tosses - heads),
    colour = "black"
  ) +
  geom_function(
    mapping = aes(linetype = "Prior"),
    fun = ~1,
    colour = "black",
    linetype = "dashed"
  ) +
  geom_freqpoly(
    aes(y = after_stat(density)),
    binwidth = 0.03,
    position = "dodge",
    linewidth = 0.75
  ) +
  scale_linetype_manual(
    values = c(
      "Prior" = "dashed",
      "True Posterior" = "solid"
    )
  ) +
  theme_bw(base_family = "serif") +
  theme(
    text = element_text(family = "serif"),
    # legend.position = "bottom",
    # legend.box="vertical",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 3)) +
  scale_x_continuous(expand = expansion(0), limits = c(0, 1)) +
  labs(
    x = "p (Probability of Heads)",
    y = "Density",
    fill = NULL,
    linetype = NULL
  )

ggsave(here("write_up", "images", "coin_MH_R.pdf"), width = 5.5, height = 2.5)
