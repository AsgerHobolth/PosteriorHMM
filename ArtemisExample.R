#######################################
##### Example of Artemis analysis #####
#######################################

# Load required packages
library(Rfast)
library(rlist)
library(dplyr)
library(ggplot2)


# Load functions
source("HybridFunctions.R")


# Define model parameters
# In this example, we use case (ii) with a=5
a <- 5

# Number of states
m <- 3

# Initial probabilities
pi <- c(0.8, 0.1, 0.1)

# Transition probabilities
Gamma <- matrix(c(0.8, 0.1, 0.1,
                  0.1, 0.8, 0.1,
                  0.1, 0.1, 0.8), nrow = 3, byrow = T)

# Emission rates
lambda <- c(20-a, 20, 20+a)

# Combine model parameters using the function model()
mod <- model(m, pi, Gamma, lambda)


# Artemis analysis
# Number of simulations
b <- 10

# Save optimal alphas, and pointwise accuracy and log joint probability for optimal alphas and all other alphas
points <- list()
acc_posts <- NULL

# This loop can be parallelized using future_lapply if it is too slow
for (i in 1:b) {
  print(sprintf("Simulation %s", i))
  
  # Simulate true hidden state sequence and observed sequence of size n based on model parameters using the function sim_data()
  dat <- sim_data(mod, n = 10^5)
  
  # Combine the observed sequence and n using the function observations()
  obs <- observations(dat$observed)
  
  # This function calculates the pointwise accuracy and log joint probability for each alpha where hybrid changes
  # and finds the alpha closest to the 45 degree angle
  res <- choose_alpha_accuracy_posterior(mod, obs, dat$true)
  
  # This is the pointwise accuracy and log joint probability of the optimal alpha
  points <- list.append(points, res$point)
  
  # This is the pointwise accuracy and log joint probability for each change point
  acc_posts <- bind_rows(acc_posts, data.frame(acc = res$acc_post[,1], log_post = res$acc_post[,2], alph = res$acc_post[,3], iteration = i) %>% mutate(best_alph = alph == res$alpha))
}


# Artemis plot for one simulation
ggplot() +
  theme_light() +
  geom_path(data = acc_posts %>% filter(iteration == 1), aes(x = acc, y = log_post), col = "#C8A1E0") +
  geom_point(aes(x = points[[1]][1], y = points[[1]][2]), size = 3, col = "#C8A1E0") +
  geom_text(aes(x = -Inf, y = -Inf), hjust = -0.05, vjust = -0.5, label = sprintf("Best alpha: %s", acc_posts %>% filter(iteration == 1) %>% filter(best_alph) %>% select(alph))) +
  labs(x = "Pointwise accuracy", y = "Log joint probability")


# Artemis plot for all 10 simulations
# Make x- and y-axis comparable between simulations
acc_posts %>%
  ungroup() %>%
  filter(alph == 0 | alph == 1) %>%
  mutate(max_x = ifelse(alph == 0, acc, 0)) %>%
  mutate(min_x = ifelse(alph == 1, acc, 0)) %>%
  mutate(max_y = ifelse(alph == 1, log_post, 0)) %>%
  mutate(min_y = ifelse(alph == 0, log_post, 0)) %>%
  select(-c(alph, acc, log_post)) %>%
  group_by(iteration) %>%
  summarize(max_x = sum(max_x),
            min_x = sum(min_x),
            max_y = sum(max_y),
            min_y = sum(min_y)) -> minmaxaxis

p <- acc_posts %>%
  left_join(minmaxaxis, by = "iteration") %>%
  group_by(iteration) %>%
  mutate(x = (acc-min_x)/(max_x-min_x)) %>%
  mutate(y = (log_post-min_y)/(max_y-min_y)) %>%
  ungroup()

# Optimal alpha based on average of 10 simulations
alphas <- p %>% filter(best_alph) %>% summarise(mean = mean(alph), sd = sd(alph))

# Plot
ggplot() +
  geom_path(data = p, aes(x = x, y = y, group = iteration), col = "#C8A1E0") +
  geom_point(data = p %>% filter(best_alph), aes(x = x, y = y), size = 3, col = "#C8A1E0") +
  theme_light() +
  theme(legend.position = 'none') +
  geom_text(aes(x = -Inf, y = -Inf), hjust = -0.05, vjust = -0.5, label = sprintf("Best alpha: %s (%s sd)", round(alphas$mean, digits = 3), round(alphas$sd, digits = 2))) +
  labs(x = "Pointwise accuracy", y = "Log joint probability")
