library(deSolve)
library(tidyverse)
library(viridis)

cust_pal <- palette.colors(4)
names(cust_pal) <- c("S", "I", "R", "E")

# Task 1 ----
set.seed(4503)
beta <- 0.4
gamma <- 1 / 4
gamma_SI_demog <- 1 / 90
omega <- 1 / 2
mu <- 0.012 # gamma_SI_demog - gamma_SI_demog^2/beta + 0.2; mu
# Aus birth rate 0.00004, rabbits 0.03
nu <- (beta * (mu - gamma_SI_demog) + gamma_SI_demog^2) / 
  (beta - gamma_SI_demog)
pop_size <- 1000 # population of Whitehorse 25000
init_inf <- 10

parameters_SIS <- c(beta = beta, gamma = gamma)
parameters_SI_demog <- c(beta = beta, gamma = gamma_SI_demog, mu = mu, nu = nu)
parameters_SEIR <- c(beta = beta, gamma = gamma, omega = omega)

state_SIS <- c(S = pop_size - init_inf, I = init_inf)
state_SI_demog <- c(S = pop_size - init_inf, I = init_inf)
state_SEIR <- c(S = pop_size - init_inf, E = 0, I = init_inf, R = 0)

ODE_fun_SIS <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dS <- gamma * I - beta * S * I / (S + I)
    dI <- beta * S * I / (S + I) - gamma * I
    # return the rate of change
    list(c(dS, dI))
  }) # end with(as.list ...
}

ODE_fun_SI_demog <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dS <- mu * (I + S) - beta * S * I / (S + I) - nu * S
    dI <- beta * S * I / (S + I) - (nu + gamma) * I
    # return the rate of change
    list(c(dS, dI))
  }) # end with(as.list ...
}

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

pand_length_SIS <- 90
pand_length_SI_demog <- 90
pand_length_SEIR <- 120
times_SIS <- seq(0, pand_length_SIS, by = 0.1)
times_SI_demog <- seq(0, pand_length_SI_demog, by = 0.1)
times_SEIR <- seq(0, pand_length_SEIR, by = 0.1)

outputs_SIS <- ode(
  y = state_SIS,
  times = times_SIS,
  func = ODE_fun_SIS,
  parms = parameters_SIS
) %>%
  as_tibble() %>%
  mutate(
    S = as.double(S),
    I = as.double(I),
    time = as.double(time)
  )
outputs_SI_demog <- ode(
  y = state_SI_demog,
  times = times_SI_demog,
  func = ODE_fun_SI_demog,
  parms = parameters_SI_demog
) %>%
  as_tibble() %>%
  mutate(
    S = as.double(S),
    I = as.double(I),
    time = as.double(time)
  )
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

long_outputs_SIS <- pivot_longer(
  outputs_SIS,
  c(S, I),
  names_to = "Compartment",
  values_to = "number"
) %>% mutate(model = "SIS")
long_outputs_SI_demog <- pivot_longer(
  outputs_SI_demog,
  c(S, I),
  names_to = "Compartment",
  values_to = "number"
) %>% mutate(model = "SI with demography")
long_outputs_SEIR <- pivot_longer(
  outputs_SEIR,
  c(S, E, I, R),
  names_to = "Compartment",
  values_to = "number"
) %>% mutate(model = "SEIR")

long_output_merged <- bind_rows(
  long_outputs_SIS, long_outputs_SI_demog, long_outputs_SEIR
) %>%
  mutate(
    model = fct_inorder(model),
    Compartment = factor(Compartment, levels = c("S", "E", "I", "R"))
  )

ODE_plot_SIS <- ggplot(
  long_outputs_SIS,
  aes(x = time, y = number, colour = Compartment)
) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_color_manual(values = cust_pal)
ODE_plot_SIS

ODE_plot_SI_demog <- ggplot(
  long_outputs_SI_demog,
  aes(x = time, y = number, colour = Compartment)
) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_color_manual(values = cust_pal)
ODE_plot_SI_demog

ODE_plot_SEIR <- ggplot(
  long_outputs_SEIR,
  aes(x = time, y = number, colour = Compartment)
) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(c(0, 0.1))) +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_color_manual(values = cust_pal)
ODE_plot_SEIR

ODE_plots <- ggplot(
  long_output_merged,
  aes(x = time, y = number, colour = Compartment)
) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(c(0, 0.01))) +
  scale_color_manual(values = cust_pal) +
  facet_grid(cols = vars(model), scales = "free_x")
ODE_plots

ggsave(
  plot = ODE_plot_SIS,
  here("../../Apps/Overleaf/M_Scimat_Thesis/images/ODE_SIS.pdf"),
  width = 3,
  height = 3
)
ggsave(
  plot = ODE_plot_SI_demog,
  here("../../Apps/Overleaf/M_Scimat_Thesis/images/ODE_SI_demog.pdf"),
  width = 3,
  height = 3
)
ggsave(
  plot = ODE_plot_SEIR,
  here("../../Apps/Overleaf/M_Scimat_Thesis/images/ODE_SEIR.pdf"),
  width = 3,
  height = 3
)
ggsave(
  plot = ODE_plots,
  here("../../Apps/Overleaf/M_Scimat_Thesis/images/ODE_plots.pdf"),
  width = 5.5,
  height = 3
)
