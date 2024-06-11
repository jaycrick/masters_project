rm(list = ls())

library(tidyverse)
library(here)
library(magrittr)

load(file = "stochastic_runs.RData")

I_per_day <- map_dbl(
  1:120,
  .f = ~ sim_SEIR %>%
    filter(t < .x) %>%
    tail(1) %$%
    I
)

list_of_days <- 1:7 * 14
obs_prevs <- I_per_day[list_of_days] / 1000

prev_data <- tibble(t = list_of_days, prevalence = obs_prevs)

obs_prevs


# Create SEIR ODE fn ------------------------------------------------------

cust_pal <- palette.colors(3)
set.seed(4503)
true_beta <- 0.4
gamma <- 1 / 4
omega <- 1 / 2
pop_size <- 1000 # population of Whitehorse 25000
init_inf <- 10

state_SEIR <- c(S = pop_size - init_inf, E = 0, I = init_inf, R = 0)
ODE_fun_SEIR <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dS <- -beta * S * I / (S + E + I + R)
    dE <- beta * S * I / (S + E + I + R) - omega * E
    dI <- omega * E - gamma * I
    dR <- gamma * I
    # return the rate of change
    list(c(dS, dE, dI, dR))
  }) # end with(as.list ...
}

pand_length_SEIR <- 120
times_SEIR <- seq(0, pand_length_SEIR, by = 0.1)

ODE_prevs <- function(beta = true_beta) {
  parameters_SEIR <- c(beta = beta, gamma = gamma, omega = omega)
  outputs_SEIR <- ode(
    y = state_SEIR,
    times = times_SEIR,
    func = ODE_fun_SEIR,
    parms = parameters_SEIR
  ) %>%
    as_tibble() %>%
    mutate(
      S = as.double(S),
      E = as.double(E),
      I = as.double(I),
      R = as.double(R),
      time = as.double(time)
    )

  return(map_dbl(
    1:120,
    .f = ~ outputs_SEIR %>%
      filter(time < .x) %>%
      tail(1) %$%
      I
  )[list_of_days] / 1000)
}

SSE_prev <- function(beta) sum((ODE_prevs(beta) - obs_prevs)^2)

NLL_prev <- function(beta) -sum(dbinom(obs_prevs*1000, 1000, ODE_prevs(beta), log = T))

LSE_beta <- optimise(SSE_prev, c(0, 1))$minimum

MLE_beta <- optimise(NLL_prev, c(0, 1))$minimum

parameters_SEIR <- c(beta = LSE_beta, gamma = gamma, omega = omega)
outputs_SEIR_LSE <- ode(
  y = state_SEIR,
  times = times_SEIR,
  func = ODE_fun_SEIR,
  parms = parameters_SEIR
) %>%
  as_tibble() %>%
  mutate(
    S = as.double(S),
    E = as.double(E),
    I = as.double(I),
    R = as.double(R),
    time = as.double(time)
  )

parameters_SEIR <- c(beta = MLE_beta, gamma = gamma, omega = omega)
outputs_SEIR_MLE <- ode(
  y = state_SEIR,
  times = times_SEIR,
  func = ODE_fun_SEIR,
  parms = parameters_SEIR
) %>%
  as_tibble() %>%
  mutate(
    S = as.double(S),
    E = as.double(E),
    I = as.double(I),
    R = as.double(R),
    time = as.double(time)
  )

MLE_LSE_plot <- ggplot(prev_data, aes(t, prevalence)) +
  geom_point(aes(color = cust_pal[1]), size = 3) +
  geom_line(aes(x = time, y = I/1000, colour = cust_pal[2]), data = outputs_SEIR_LSE, linewidth = 0.75) +
  geom_line(aes(x = time, y = I/1000, colour = cust_pal[3]), data = outputs_SEIR_MLE, linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(0), limits = c(0, 120)) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 0.07)) +
  xlab("Time in Days Since Start of Epidemic") +
  ylab("Disease Prevalence (Proportion)") +
  scale_color_manual(name = NULL, values = cust_pal, labels = c("Observed Data", "LSE ODE solution", "MLE ODE solution"))

ggsave("write_up/images/MLE_SSE_SEIR.pdf", MLE_LSE_plot, width = 5.5, height = 3)
