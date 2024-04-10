library(tidyverse)
library(latex2exp)

x <- 1:3
y <- c(2, 4, 4)

mod <- lm(y ~ x)

pois_NLL <- function(params) {
  a <- params[[1]]
  b <- params[[2]]
  NLL <- -dpois(2, a + b, log = T) - dpois(4, a + 3 * b, log = T) - dpois(4, a + 4 * b, log = T)
  return(NLL)
}

pois_pars <- optim(c(1, 5 / 4), pois_NLL)$par

ggplot(mapping = aes(x = x, y = y)) +
  geom_point(aes(fill = "Observations"), size = 2) +
  theme_bw() +
  geom_abline(
    aes(
      slope = mod$coefficients[[2]],
      intercept = mod$coefficients[[1]],
      colour = "Least Squares Estimate"
    )
  ) +
  geom_abline(
    aes(
      slope = pois_pars[[2]],
      intercept = pois_pars[[1]],
      colour = "Poisson Maximum\nLikelihood Estimate"
    )
  ) +
  theme(
    text = element_text(family = "serif"),
    # legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(expand = expansion(0), limits = c(0, 3.5)) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 4.5))

ggsave("C:/Users/jckricket/Dropbox/Apps/Overleaf/M_Scimat_Thesis/images/LS_example.pdf", width = 6, height = 2.5)
