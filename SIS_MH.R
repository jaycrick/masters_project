library(deSolve)
library(latex2exp)
library(tidyverse)

# Task 1 ----
set.seed(4503)
SIS_beta <- 0.4
SIS_gamma <- 1 / 4 # from rgamma(1, 2, 6)
pop_size <- 1000 # population of Whitehorse 25000

parameters <- c(beta = SIS_beta, gamma = SIS_gamma)
state <- c(s = (1 - 1 / pop_size), i = 1 / pop_size)

SIS_function <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    ds <- gamma * i - beta * s * i
    di <- beta * s * i - gamma * i
    # return the rate of change
    list(c(ds, di))
  }) # end with(as.list ...
}

pand_length <- 60

times <- seq(0, pand_length, by = 0.1)

outputs <- ode(
  y = state,
  times = times,
  func = SIS_function,
  parms = parameters
) %>%
  as_tibble() %>%
  mutate(
    s = as.double(s),
    i = as.double(i),
    time = as.double(time)
  )

long_outputs <- pivot_longer(
  outputs, c(s, i),
  names_to = "state",
  values_to = "proportions"
)

I_t <- function(
    t,
    b = SIS_beta,
    g = SIS_gamma,
    population = pop_size,
    initial_infections = 1) {
  if (b == g) {
    return(initial_infections)
  }
  i_inf <- (1 - g / b) * population
  chi <- b - g
  V <- i_inf / initial_infections - 1
  i_inf / (1 + V * exp(-chi * t))
}

SIS_plot <- ggplot(
  long_outputs,
  aes(x = time, y = proportions, colour = state)
) +
  geom_line(linewidth = 0.75) +
  theme_bw() +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since first infection") +
  ylab("Proportion in each compartment") +
  geom_function(aes(linetype = "dashed"), fun = ~ I_t(.x) / pop_size)
SIS_plot

sample_day <- 30

sample_day_prop_i <- long_outputs %>%
  filter(time == sample_day & state == "i") %>%
  select(state, proportions)

set.seed(4503)
case_estimate <- rpois(
  1,
  SIS_beta * I_t(sample_day) * (pop_size - I_t(sample_day)) / pop_size
)

binom_rho_accept <- function(
    x,
    y,
    cases = case_estimate,
    population = pop_size) {
  log_likihood_ratio <- dbinom(
    cases,
    (population - round(I_t(sample_day, b = y), 0)),
    SIS_beta * I_t(sample_day, b = y) / population,
    log = T
  ) + dgamma(y, 2, 6, log = T) -
    dbinom(
      cases,
      (population - round(I_t(sample_day, b = x), 0)),
      SIS_beta * I_t(sample_day, b = x) / population,
      log = T
    ) -
    dgamma(x, 2, 6, log = T)
  likelihood_ratio <- exp(log_likihood_ratio)
  min(likelihood_ratio, 1)
}

no_samples <- 30000
init_sample <- 0.5

binom_MH_samples <- c(init_sample)

for (i in 2:no_samples) {
  last_sample <- binom_MH_samples[i - 1]
  new_guess <- rnorm(1, mean = last_sample, sd = 1 / 100)
  accept_prob <- binom_rho_accept(last_sample, new_guess)
  binom_MH_samples[i] <- sample(c(last_sample, new_guess), 1, prob = c(1 - accept_prob, accept_prob))
}

binom_MH_tbl <- tibble(Samples = binom_MH_samples[1001:no_samples])

# gamma_dens <- dgamma(1/4 + (1:999) / 100000, 2, 6, log = T)
#
# infecteds_list <- sapply(0.225 + (1:999) / 100000, function(x) I_t(sample_day, b = x))
#
# binom_infected_likelihoods <- sapply(infecteds_list, function(x) dbinom(case_estimate, pop_size - round(x, 0), SIS_beta * x / pop_size, log = T))
#
# binom_posteriors <- gamma_dens + binom_infected_likelihoods
# true_post_binom <- tibble(x = 0.225 + (1:999) / 100000, y = exp(binom_posteriors) / (sum(exp(binom_posteriors)) / 100000))

ggplot(binom_MH_tbl, mapping = aes(x = Samples, fill = "Sampled posterior")) +
  geom_histogram(aes(y = after_stat(density)), bins = 40) +
  # geom_line(data = true_post_binom, mapping = aes(x = x, y = y, colour = "True posterior")) +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = expansion(0)) +
  scale_x_continuous(expand = expansion(0)) +
  labs(title = "Samples from beta", x = "beta", y = "Density", colour = NULL, fill = NULL) +
  theme_bw()

pois_rho_accept <- function(x, y, cases = case_estimate, population = pop_size) {
  log_likihood_ratio <- dpois(
    cases,
    SIS_beta *
      I_t(sample_day, b = y) *
      (population - I_t(sample_day, b = y)) /
      population,
    log = T
  ) +
    dgamma(y, 2, 6, log = T) -
    dpois(
      cases,
      SIS_beta *
        I_t(sample_day, b = x) *
        (population - I_t(sample_day, b = x)) /
        population,
      log = T
    ) -
    dgamma(x, 2, 6, log = T)
  likelihood_ratio <- exp(log_likihood_ratio)
  min(likelihood_ratio, 1)
}

no_samples <- 30000
init_sample <- 0.23

pois_MH_samples <- c(init_sample)

for (i in 2:no_samples) {
  last_sample <- pois_MH_samples[i - 1]
  new_guess <- rnorm(1, mean = last_sample, sd = 1 / 1000)
  accept_prob <- pois_rho_accept(last_sample, new_guess)
  pois_MH_samples[i] <- sample(
    c(last_sample, new_guess),
    1,
    prob = c(1 - accept_prob, accept_prob)
  )
}

pois_MH_tbl <- tibble(Samples = pois_MH_samples[1001:no_samples])

# pois_infected_likelihoods <- sapply(infecteds_list, function(x) dpois(case_estimate, SIS_beta * x * (pop_size - x) / pop_size, log = T))
#
# pois_posteriors <- gamma_dens + pois_infected_likelihoods
# true_post_pois <- tibble(x = 0.225 + (1:999) / 100000, y = exp(pois_posteriors) / (sum(exp(pois_posteriors)) / 100000))

ggplot(pois_MH_tbl, mapping = aes(x = Samples, fill = "Sampled posterior")) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, position = "dodge") +
  # geom_line(data = true_post_pois, mapping = aes(x = x, y = y, colour = "Poisson posterior")) +
  theme(
    text = element_text(family = "serif"),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_y_continuous(expand = expansion(0)) +
  scale_x_continuous(expand = expansion(0)) +
  labs(title = "Samples from beta", x = "beta", y = "Density", colour = NULL, fill = NULL) +
  theme_bw()

both_MH_tbl <- tibble(
  `Binomial Likelihood` = binom_MH_samples[1001:no_samples],
  `Poisson Likelihood` = pois_MH_samples[1001:no_samples]
) %>%
  pivot_longer(
    everything(),
    names_to = "Sampled Posterior"
  )

ggplot(both_MH_tbl, mapping = aes(x = value)) +
  geom_histogram(aes(y = after_stat(density), fill = `Sampled Posterior`), bins = 40, colour = "black") +
  # geom_line(data = true_post_pois, mapping = aes(x = x, y = y, linetype = "Poisson")) +
  # geom_line(data = true_post_binom, mapping = aes(x = x, y = y, linetype = "Binomial")) +
  theme_bw() +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(0, 1)) +
  theme(
    text = element_text(family = "serif"),
    # legend.position = "bottom",
    # legend.box="vertical",
    # plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    x = TeX(r"($\beta$)"),
    y = "Density",
    colour = NULL,
    # linetype = "True Posterior"
  ) +
  guides(fill = "none") +
  facet_wrap(vars(`Sampled Posterior`), ncol = 2)

ggsave("write_up/images/SIS_beta_pred.pdf", width = 5.5, height = 3)
