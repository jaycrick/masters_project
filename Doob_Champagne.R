library(tidyverse)
library(here)


pv_champ_alpha <- 0.4 # prop of effective care
pv_champ_beta <- 0.4 # prop of radical cure
pv_champ_gamma_L <- 1 / 223 # liver stage clearance rate
pv_champ_delta <- 0.05 # prop of imported cases
pv_champ_lambda <- 0.04 # transmission rate
pv_champ_f <- 1 / 72 # relapse frequency
pv_champ_r <- 1 / 60 # blood stage clearance rate


champagne_stochastic <- function(alpha, beta, gamma_L, lambda, f, r, N = 1000, I_L = 10, I_0 = 0, S_L = 0, delta = 0, end_time = 1000) {
  t <- 0
  S_0 <- N - I_L - I_0 - S_L
  list_of_outcomes <- list(c(t = 0, S_0 = S_0, S_L = S_L, I_0 = I_0, I_L = I_L))
  while (t < end_time) {
    if (S_0 == N) {
      break
    }
    S_0_to_I_L <- (1 - alpha) * lambda * (I_L + I_0) / N * S_0 
    S_0_to_S_L <- alpha * (1 - beta) * lambda * (I_0 + I_L) / N * S_0
    I_0_to_S_0 <- r * I_0 / N
    I_0_to_I_L <- lambda * (I_L + I_0) / N * I_0
    I_L_to_I_0 <- gamma_L * I_L
    I_L_to_S_L <- r * I_L
    S_L_to_S_0 <- (gamma_L + (f + lambda * (I_0 + I_L) / N) * alpha * beta) * S_L
    S_L_to_I_L <- (f + lambda * (I_0 + I_L) / N) * (1 - alpha) * S_L
    total_rate <- S_0_to_I_L + S_0_to_S_L +
      I_0_to_S_0 + I_0_to_I_L +
      I_L_to_I_0 + I_L_to_S_L +
      S_L_to_S_0 + S_L_to_I_L
    # total_rate <- (lambda * (I_0 + I_L) + delta) * alpha * beta * S_L + alpha*beta*f*S_L + gamma_L*S_L + r*I_0 + #entering S_0
    # gamma_L*I_L + #Entering I_0
    # alpha*(1-beta)*(lambda*(I_0+I_L) + delta)*S_0 + r*I_L + #entering S_L
    # (1 - alpha)*(lambda*(I_L + I_0) + delta)*(S_0 + S_L) + (lambda*(I_L + I_0) + delta)*I_0 + (1 - alpha)*f*S_L #entering I_L
    t <- t + rexp(1, total_rate)
    new_stages <- sample(
      list(
        c(t = t, S_0 = S_0 - 1, S_L = S_L, I_0 = I_0, I_L = I_L + 1), c(t = t, S_0 = S_0 - 1, S_L = S_L + 1, I_0 = I_0, I_L = I_L),
        c(t = t, S_0 = S_0 + 1, S_L = S_L, I_0 = I_0 - 1, I_L = I_L), c(t = t, S_0 = S_0, S_L = S_L, I_0 = I_0 - 1, I_L = I_L + 1),
        c(t = t, S_0 = S_0, S_L = S_L, I_0 = I_0 + 1, I_L = I_L - 1), c(t = t, S_0 = S_0, S_L = S_L + 1, I_0 = I_0, I_L = I_L - 1),
        c(t = t, S_0 = S_0 + 1, S_L = S_L - 1, I_0 = I_0, I_L = I_L), c(t = t, S_0 = S_0, S_L = S_L - 1, I_0 = I_0, I_L = I_L + 1)
      ),
      1,
      prob = c(
        S_0_to_I_L / total_rate, S_0_to_S_L / total_rate,
        I_0_to_S_0 / total_rate, I_0_to_I_L / total_rate,
        I_L_to_I_0 / total_rate, I_L_to_S_L / total_rate,
        S_L_to_S_0 / total_rate, S_L_to_I_L / total_rate
      )
    )
    list_of_outcomes <- append(list_of_outcomes, new_stages)

    S_0 <- new_stages[[1]][["S_0"]]
    I_0 <- new_stages[[1]][["I_0"]]
    I_L <- new_stages[[1]][["I_L"]]
    S_L <- new_stages[[1]][["S_L"]]
  }
  outcome_tbl <- tibble(list_of_outcomes) %>%
    unnest_wider(list_of_outcomes, names_repair = ~ c("t", "S_0", "I_0", "I_L", "S_L")) %>%
    suppressMessages()
  return(outcome_tbl)
}

set.seed(590154)

champ_samp <- champagne_stochastic(pv_champ_alpha, pv_champ_beta, pv_champ_gamma_L, pv_champ_lambda, pv_champ_f, pv_champ_r) %>%
  pivot_longer(-t)

ggplot(champ_samp, aes(x = t, y = value, colour = name)) +
  geom_line() +
  theme_classic() +
  scale_x_continuous(expand = expansion(0)) +
  scale_y_continuous(expand = expansion(0))

ggsave("C:/Users/jckricket/Dropbox/Apps/Overleaf/M_Scimat_Thesis/images/doob_champagne.pdf", width = 6, height = 5, dpi = 600)
