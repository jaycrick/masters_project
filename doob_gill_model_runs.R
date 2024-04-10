library(tidyverse)
library(here)

cust_pal <- palette.colors(4)
names(cust_pal) <- c("S", "I", "R", "E")

# Task 1 SIR ----

set.seed(4503)
beta <- 0.4
gamma <- 1 / 4
gamma_SI_demog <- 1 / 90
omega <- 1 / 2
mu <- 0.012 # gamma_SI_demog - gamma_SI_demog^2/beta + 0.001
# Aus birth rate 0.00004
nu <- (beta * (mu - gamma_SI_demog) + gamma_SI_demog^2) / 
  (beta - gamma_SI_demog)
pop_size <- 1000 # population of Whitehorse 25000
init_inf <- 10

S_0 <- pop_size - init_inf
E_0 <- 0
I_0 <- init_inf
R_0 <- 0

sim_SIS <- tibble(t = 0, S = pop_size - init_inf, I = init_inf)
sim_SI_demog <- tibble(t = 0, S = pop_size - init_inf, I = init_inf)
sim_SEIR <- tibble(
  t = 0, S = pop_size - init_inf, E = 0, I = init_inf, R = 0
)

S_current <- S_0
I_current <- I_0
t_current <- 0
while (I_current > 0 & t_current < 90) {
  inf_rate <- beta * S_current * I_current / (pop_size - 1)
  rec_rate <- gamma * I_current
  tot_rate <- inf_rate + rec_rate
  t_current <- t_current + rexp(1, tot_rate)

  if (rbinom(1, 1, inf_rate / tot_rate)) {
    new_state <- tibble(
      "t" = c(t_current, t_current),
      "S" = c(S_current, S_current - 1),
      "I" = c(I_current, I_current + 1)
    )
  } else {
    new_state <- tibble(
      "t" = c(t_current, t_current),
      "S" = c(S_current, S_current + 1),
      "I" = c(I_current, I_current - 1)
    )
  }
  sim_SIS <- bind_rows(sim_SIS, new_state)
  S_current <- unlist(new_state["S"])[[2]]
  I_current <- unlist(new_state["I"])[[2]]
}

S_current <- S_0
I_current <- I_0
t_current <- 0
while (I_current > 0 & t_current < 45) {
  birth_rate <- mu * (S_current + I_current)
  death_rate <- nu * S_current
  inf_rate <- beta * S_current * I_current / (I_current + S_current)
  rec_rate <- (gamma_SI_demog + nu) * I_current
  tot_rate <- birth_rate + death_rate + inf_rate + rec_rate
  t_current <- t_current + rexp(1, tot_rate)

  event <- sample(
    c(1, 2, 3, 4),
    1,
    prob = c(birth_rate, death_rate, inf_rate, rec_rate) / tot_rate,
    replace = T
  )
  new_state <- case_when(
    event == 1 ~
      tibble(
        "t" = c(t_current, t_current),
        "S" = c(S_current, S_current + 1),
        "I" = c(I_current, I_current)
      ),
    event == 2 ~
      tibble(
        "t" = c(t_current, t_current),
        "S" = c(S_current, S_current - 1),
        "I" = c(I_current, I_current)
      ),
    event == 3 ~
      tibble(
        "t" = c(t_current, t_current),
        "S" = c(S_current, S_current - 1),
        "I" = c(I_current, I_current + 1)
      ),
    event == 4 ~
      tibble(
        "t" = c(t_current, t_current),
        "S" = c(S_current, S_current),
        "I" = c(I_current, I_current - 1)
      )
  )
  sim_SI_demog <- bind_rows(sim_SI_demog, new_state)
  S_current <- unlist(new_state["S"])[[2]]
  I_current <- unlist(new_state["I"])[[2]]
}

S_current <- S_0
E_current <- E_0
I_current <- I_0
R_current <- R_0
t_current <- 0
while (I_current + E_current > 0 & t_current < 120) {
  inf_rate <- beta * S_current * I_current / pop_size
  lat_rate <- omega * E_current
  rec_rate <- gamma * I_current
  tot_rate <- inf_rate + lat_rate + rec_rate
  t_current <- t_current + rexp(1, tot_rate)

  event <- sample(
    c(1, 2, 3),
    1,
    prob = c(inf_rate, lat_rate, rec_rate) / tot_rate,
    replace = T
  )
  new_state <- case_when(
    event == 1 ~
      tibble(
        "t" = c(t_current, t_current),
        "S" = c(S_current, S_current - 1),
        "E" = c(E_current, E_current + 1),
        "I" = c(I_current, I_current),
        "R" = c(R_current, R_current)
      ),
    event == 2 ~
      tibble(
        "t" = c(t_current, t_current),
        "S" = c(S_current, S_current),
        "E" = c(E_current, E_current - 1),
        "I" = c(I_current, I_current + 1),
        "R" = c(R_current, R_current)
      ),
    event == 3 ~
      tibble(
        "t" = c(t_current, t_current),
        "S" = c(S_current, S_current),
        "E" = c(E_current, E_current),
        "I" = c(I_current, I_current - 1),
        "R" = c(R_current, R_current + 1)
      )
  )
  sim_SEIR <- bind_rows(sim_SEIR, new_state)
  S_current <- unlist(new_state["S"])[[2]]
  E_current <- unlist(new_state["E"])[[2]]
  I_current <- unlist(new_state["I"])[[2]]
  R_current <- unlist(new_state["R"])[[2]]
}

long_outputs_SIS <- pivot_longer(
  sim_SIS,
  c(S, I),
  names_to = "Compartment",
  values_to = "number"
) %>% mutate(model = "SIS")

long_outputs_SI_demog <- pivot_longer(
  sim_SI_demog,
  c(S, I),
  names_to = "Compartment",
  values_to = "number"
) %>% mutate(model = "SI with demography")

long_outputs_SEIR <- pivot_longer(
  sim_SEIR,
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

plot_stochastic_SIS <- ggplot(
  long_outputs_SIS,
  aes(x = t, y = number, colour = Compartment)
) +
  geom_path(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_color_manual(values = cust_pal)

plot_stochastic_SI_demog <- ggplot(
  long_outputs_SI_demog,
  aes(x = t, y = number, colour = Compartment)
) +
  geom_path(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_color_manual(values = cust_pal)

plot_stochastic_SEIR <- ggplot(
  long_outputs_SEIR,
  aes(x = t, y = number, colour = Compartment)
) +
  geom_path(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(0)) +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_color_manual(values = cust_pal)

plots_stochastic_all <- ggplot(
  long_output_merged,
  aes(x = t, y = number, colour = Compartment)
) +
  geom_line(linewidth = 0.75) +
  theme_bw(base_family = "serif") +
  xlab("Days since pandemic start") +
  ylab("Number of people") +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(
    expand = expansion(c(0, 0.01)), 
    breaks = c(0, 250, 500, 750, 1000, 1250)
  ) +
  scale_color_manual(values = cust_pal) +
  facet_grid(cols = vars(model), scales = "free_x")
plots_stochastic_all

ggsave(
  plot = plot_stochastic_SIS,
  "../../Apps/Overleaf/M_Scimat_Thesis/images/doob_SIS.pdf",
  width = 3,
  height = 3
)
ggsave(
  plot = plot_stochastic_SI_demog,
  "../../Apps/Overleaf/M_Scimat_Thesis/images/doob_SI_demog.pdf",
  width = 3,
  height = 3
)
ggsave(
  plot = plot_stochastic_SEIR,
  "../../Apps/Overleaf/M_Scimat_Thesis/images/doob_SEIR.pdf",
  width = 3,
  height = 3
)
ggsave(
  plot = plots_stochastic_all,
  here("../../Apps/Overleaf/M_Scimat_Thesis/images/doob_plots.pdf"),
  width = 5.5,
  height = 3
)
