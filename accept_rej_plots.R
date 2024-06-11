rm(list = ls())

library(tidyverse)

samples = 500

set.seed(31097409)

my_tbl = tibble(x = runif(samples, 0, 2), y =  runif(samples, 0, 1)) %>%
  mutate(
    Colour = if_else(y < (x - 1)^2, "Accept", "Reject")
  )

table(my_tbl$Colour)

accept_rej_plot <- ggplot(my_tbl, aes(x, y, color = Colour)) +
  geom_point(size = 2) +
  theme_bw(base_family = "serif") +
  scale_x_continuous(expand = expansion(0), limits = c(0, 2)) +
  scale_y_continuous(expand = expansion(0), limits = c(0, 1)) +
  geom_function(aes(colour = "g(x) = (x - 1)^2"), fun = ~(.x - 1)^2, linewidth = 0.75) +
  scale_color_manual(name = NULL, values = c("3", "1", "2"), labels = c("Accepted Samples", "Unnormalised\nProbability Density", "Rejected Samples"))+
  xlab("X*") +
  ylab("U")

accept_rej_plot

ggsave("write_up/images/accept_reject.pdf", accept_rej_plot, width = 5.5, height = 3)
