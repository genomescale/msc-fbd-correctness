library(scales)
library(ggplot2)
library(cowplot)

estimated_marginals = read.csv("estimated-marginals.csv")
estimated_marginals$config = as.factor(estimated_marginals$config)
estimated_marginals$variable = as.factor(estimated_marginals$variable)

topology0_marginals = subset(estimated_marginals, topology_code == 0)

true_marginals = read.csv("true-marginals.csv")
true_marginals$variable = as.factor(true_marginals$variable)

levels(topology0_marginals$variable) = c(expression("t"["or"]), expression("x"["1"]), expression("x"["2"]))
levels(true_marginals$variable) = c(expression("t"["or"]), expression("x"["1"]), expression("x"["2"]))
levels(topology0_marginals$config) = c("Coordinated", "Full", "MSC", "NodeReheight2", "UpDown", "SA")

ggplot(topology0_marginals, aes(x = lower_bound, y = marginal_probability * 10)) +
  geom_step(alpha = 0.1, aes(group = rep_i), color = "green") +
  geom_step(data = true_marginals) +
  facet_grid(config ~ variable, labeller = label_parsed) +
  scale_x_continuous(limits = c(0, 5)) +
  labs(x = "Time", y = "Marginal probability") +
  theme_minimal_grid() +
  theme(legend.position = "None")

ggsave("marginal-probability.png", units = "in", width = 8, height = 10)

