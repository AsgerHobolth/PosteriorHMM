model <- function(m, pi, Gamma, lambda) {
  # A function to collect model parameters for a Poisson HMM.
  #
  # Arguments:
  #   m: number of states
  #   n: number of observations
  #   observed: observed sequence (vector of length n)
  #   pi: initial distribution (vector of length m)
  #   Gamma: transition probabilities (mxm matrix)
  #   lambda: Poisson rates (vector of length m)
  #
  # Returns:
  #   A list of model parameters.
  
  result <- list(m = m, pi = pi, Gamma = Gamma, lambda = lambda)
  return(result)
}

observations <- function(observed, n = length(observed)) {
  # A function to collect the observed sequence and the number of observations n.
  #
  # Arguments:
  #   observed: observed sequence (vector of length n)
  #   n: number of observations
  #
  # Returns:
  #   A list of the observed sequence and the number of observations n.
  
  result <- list(observed = observed, n = n)
  return(result)
}

log_forward <- function(mod, obs) {
  # A function to calculate the forward probabilities and log likelihood of data.
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #
  # Returns:
  #   alpha: forward probabilities (mxn matrix)
  #   log_likelihood: log likelihood of data
  
  # create empty matrix of forward probabilities
  lalpha <- matrix(0, nrow = mod$m, ncol = obs$n)
  
  # first column in lalpha
  foo <- mod$pi * dpois(obs$observed[1], mod$lambda)
  sumfoo <- sum(foo)
  # to calculate likelihood
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[,1] <- lscale + log(foo)
  
  # rest of the columns
  for (i in 2:obs$n) {
    foo <- foo %*% mod$Gamma * dpois(obs$observed[i], mod$lambda)
    sumfoo <- sum(foo)
    # to calculate likelihood
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[,i] <- lscale + log(foo)
  }
  
  result <- list(alpha = lalpha, log_likelihood = lscale)
  return(result)
}

log_backward <- function(mod, obs) {
  # A function to calculate the backward probabilities.
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #
  # Returns:
  #   beta: backward probabilities (mxn matrix)
  
  # create empty matrix of backward probabilities
  lbeta <- matrix(0, nrow = mod$m, ncol = obs$n)
  
  lbeta[,obs$n] <- rep(0, mod$m)
  foo <- rep(1 / mod$m, mod$m)
  lscale <- log(mod$m)
  
  for (i in (obs$n-1):1) {
    foo <- mod$Gamma %*% (dpois(obs$observed[i+1], mod$lambda) * foo)
    lbeta[,i] <- log(foo) + lscale
    sumfoo <- sum(foo)
    foo <- foo / sumfoo
    lscale <- lscale + log(sumfoo)
  }
  
  return(lbeta)
}

forward_backward_llk <- function(mod, obs) {
  # A function to collect the forward and backward probabilities, and the log likelihood of data.
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #
  # Returns:
  #   A list of the forward probabilities, backward probabilities, and log likelihood of data.
  
  forw <- log_forward(mod, obs)
  back <- log_backward(mod, obs)
  
  result <- list(alpha = forw$alpha, beta = back, log_likelihood = forw$log_likelihood)
  return(result)
}

state_probs <- function(mod, n, fb_llk) {
  # A function to calculate the state probabilities P(Y_t = j | X^(n) = x^(n)).
  #
  # Arguments:
  #   mod: model parameters
  #   n: number of observations
  #   fb_llk: forward and backward probabilities, and log likelihood of data
  #
  # Returns:
  #   probs: state probabilities (mxn matrix)
  
  probs <- matrix(0, nrow = mod$m, ncol = n)
  
  for (i in 1:n) {
    probs[,i] <- exp(fb_llk$alpha[,i] + fb_llk$beta[,i] - fb_llk$log_likelihood)
  }
  
  return(probs)
}

log_joint_prob <- function(mod, obs, path) {
  # A function to calculate the log joint probability P(Y^(n) = y^(n), X^(n) = x^(n)).
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #   path: (estimated) hidden state sequence
  #
  # Returns:
  #   log_p: log joint probability
  
  log_p <- log(mod$pi[path[1]]) + log(dpois(obs$observed[1], mod$lambda[path[1]]))
  
  for (i in 2:obs$n) {
    log_p <- log_p + log(mod$Gamma[path[i-1],path[i]]) + log(dpois(obs$observed[i], mod$lambda[path[i]]))
  }
  
  return(log_p)
}

accuracy <- function(true, path) {
  # A function to calculate the pointwise accuracy of a given path.
  #
  # Arguments:
  #   true: true hidden state sequence
  #   path: (estimated) hidden state sequence
  #
  # Returns:
  #   Pointwise accuracy
  
  return(sum(path == true) / length(true)) 
}

posterior_decoding <- function(mod, n, fb_llk) {
  # A function to calculate the Posterior decoding path.
  #
  # Arguments:
  #   mod: model parameters
  #   n: number of observations
  #   fb_llk: forward and backward probabilities, and log likelihood of data
  #
  # Returns:
  #   path: Posterior decoding path
  
  probs <- state_probs(mod, n, fb_llk)
  path <- apply(probs, 2, which.max)
  
  return(path)
}

viterbi <- function(mod, obs) {
  # A function to calculate the Viterbi path.
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #
  # Returns:
  #   path: Viterbi path
  
  xi <- matrix(0, nrow = mod$m, ncol = obs$n)
  
  foo <- mod$pi * dpois(obs$observed[1], mod$lambda)
  xi[,1] <- foo / sum(foo)
  
  for (i in 2:obs$n) {
    foo <- apply(xi[,i-1] * mod$Gamma, 2, max) * dpois(obs$observed[i], mod$lambda)
    xi[,i] <- foo / sum(foo)
  }
  
  path <- rep(0, obs$n)
  
  path[obs$n] <- which.max(xi[,obs$n])
  
  for (i in (obs$n-1):1) {
    path[i] <- which.max(mod$Gamma[,path[i+1]] * xi[,i])
  }
  
  return(path)
}

hybrid <- function(mod, obs, fb_llk, alpha) {
  # A function to calculate the hybrid path for a given alpha.
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #   fb_llk: forward and backward probabilities, and log likelihood of data
  #   alpha: tuning parameter alpha
  #
  # Returns:
  #   path: hybrid path
  
  # P(Y_t = j | X^(n) = x^(n))
  phi <- state_probs(mod, obs$n, fb_llk)
  
  delta <- matrix(0, nrow = mod$m, ncol = obs$n)
  psi <- matrix(0, nrow = mod$m, ncol = obs$n-1)
  
  # create matrix of emission probabilities
  dpois_values <- matrix(dpois(rep(obs$observed, each = mod$m), rep(mod$lambda, times = obs$n)), nrow = mod$m)
  
  # calculate first column of delta
  delta[,1] <- ifelse(mod$pi == 0, (1-alpha) * log(phi[,1]), alpha * log(mod$pi * dpois_values[,1]) + (1-alpha) * log(phi[,1]))
  
  # calculate delta and psi
  for (i in 2:obs$n) {
    prob <- eachrow(y = delta[,i-1], x = alpha * log(t(mod$Gamma) * dpois_values[,i]), "+")
    delta[,i] <- rowMaxs(prob, TRUE) + (1-alpha) * log(phi[,i])
    psi[,i-1] <- rowMaxs(prob, FALSE)
  }
  
  # backtracking
  yhat <- rep(0, obs$n)
  yhat[obs$n] <- which.max(delta[,obs$n])
  for (i in (obs$n-1):1) {
    yhat[i] <- psi[yhat[i+1],i]
  }
  
  return(yhat)
}

hybrid_curry <- function(mod, obs) {
  # A function to create a hybrid function with a fixed model and observed sequence (only depends on alpha).
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #
  # Returns:
  #   A hybrid function with a fixed model and observed sequence.
  
  fb_llk <- forward_backward_llk(mod, obs)
  
  return(
    function(alpha) {
      return(hybrid(mod, obs, fb_llk, alpha))
    }
  )
}

hybrid_change <- function(hc, max_depth = 8, from = 0, to = 1, from_hyb = NULL) {
  # Function to find the change points in the hybrid path for a given hybrid curry function, using a recursive divide and conquer.
  #
  # Arguments:
  #   hc: hybrid curry function
  #   max_depth: number of times we divide the interval
  #
  # Returns:
  #   A list of change points in the hybrid path. The list is empty if no change points are found.
  
  # end of recursion
  if (max_depth < 1) {
    return(c())
  }
  
  # from hybrid path
  if (is.null(from_hyb)) {
    from_hyb <- hc(from)
  }
  
  # initial values
  change_points <- c()
  
  # middle point
  middle <- (from + to) / 2
  
  # calculate hybrid path
  hyb <- hc(middle)
  
  # is all entries from hyb equal to from_hyb
  is_equal <- all(hyb == from_hyb)
  
  # find change points
  if (!is_equal) {
    change_points <- c(change_points, list(list(alpha = middle, hyb = hyb)))
    change_points <- c(change_points, hybrid_change(hc, max_depth - 1, from, middle, from_hyb))
  }
  change_points <- c(change_points, hybrid_change(hc, max_depth - 1, middle, to, hyb))
  
  # return change points
  if (length(change_points) == 0) {
    return(c())
  }
  
  # sort change points
  change_points <- change_points[sort.list(sapply(change_points, function(x) x$alpha))]
  
  # remove duplicates
  if (length(change_points) > 1) {
    new_change_points <- c(list(change_points[[1]]))
    for (i in 2:length(change_points)) {
      if (!all(change_points[[i]]$hyb == change_points[[i-1]]$hyb)) {
        new_change_points <- c(new_change_points, list(change_points[[i]]))
      }
    }
    change_points <- new_change_points
  }
  
  return(change_points)
}

sim_data <- function(mod, n, num_states = mod$m) {
  # Function to simulate a true hidden state sequence and an observed sequence from a given model.
  #
  # Arguments:
  #   mod: model parameters
  #   n: number of observations
  #   num_states: number of states
  #
  # Returns:
  #   A list of the true hidden state sequence and the observed sequence.
  
  # generate true hidden sequence
  true <- rep(0, n)
  true[1] <- sample(1:num_states, size = 1, prob = mod$pi)
  for (i in 2:n) {
    true[i] <- sample(1:num_states, size = 1, prob = mod$Gamma[true[i-1],])
  }
  
  # generate observations
  simulated <- replicate(n, 0)
  for (i in 1:n) {
    simulated[i] <- rpois(1, mod$lambda[true[i]])
  }
  
  result <- list(true = true, observed = simulated)
  return(result)
}

accuracy_posterior <- function(mod, obs, true) {
  # A function to calculate the pointwise accuracy and posterior for each change point.
  #
  # Arguments:
  #   mod: model parameters
  #   obs: observed sequence and number of observations
  #   true: true hidden state sequence
  #
  # Returns:
  #   A matrix of pointwise accuracy and posterior for each change point.
  
  # find change points
  changes <- hybrid_change(hybrid_curry(mod, obs))
  
  posterior <- posterior_decoding(mod, obs$n, forward_backward_llk(mod, obs))
  acc_post <- c(accuracy(true, posterior), log_joint_prob(mod, obs, posterior), 0)
  
  if (length(changes) > 0) {
    for (i in 1:length(changes)) {
      acc_post <- rbind(acc_post,
                        c(accuracy(true, changes[[i]]$hyb), log_joint_prob(mod, obs, changes[[i]]$hyb), changes[[i]]$alpha))
    }
  }
  
  vit <- viterbi(mod, obs)
  acc_post <- rbind(acc_post,
                    c(accuracy(true, vit), log_joint_prob(mod, obs, vit), 1))
  
  colnames(acc_post) <- c("accuracy", "posterior", "alpha")
  
  return(acc_post)
}

distance_to_point_angle <- function(start_point = c(0, 0), angle, point) {
  distance <- abs(cos(angle) * (start_point[2] - point[2]) - sin(angle) * (start_point[1] - point[1]))
  
  return(list(distance = distance, point = point))
}


choose_alpha_accuracy_posterior <- function(mod = NULL, obs = NULL, true = NULL, angle = 45, acc_post = NULL) {
  if (is.null(acc_post)) {
    acc_post <- accuracy_posterior(mod, obs, true)
  }
  
  # transform axes to [0,1]
  acc_01 <- (acc_post[,1] - min(acc_post[,1])) / (max(acc_post[,1]) - min(acc_post[,1]))
  post_01 <- (acc_post[,2] - min(acc_post[,2])) / (max(acc_post[,2]) - min(acc_post[,2]))
  
  # transform angle to radians, since slope of line is tan(angle) (angle is in radians)
  angle <- angle * base::pi / 180
  
  # loop over all change points and find point with minimum distance to line
  # that starts in (0,0) with angle
  distances <- c()
  ps <- list()
  for (i in 2:(length(acc_post[,1])-1)) {
    s <- distance_to_point_angle(angle = angle, point = c(acc_01[i], post_01[i]))
    
    distances <- c(distances, s$distance)
    ps <- list.append(ps, s$point)
  }
  
  min_distance_index <- which.min(distances)
  optimal_point_01 <- ps[[min_distance_index]]
  r <- acc_post[2:(length(acc_post[,1])-1),]
  optimal_point <- r[min_distance_index,1:2]
  optimal_alpha <- acc_post[min_distance_index,3]
  
  res <- list(point = optimal_point, point_01 = optimal_point_01, alpha = optimal_alpha, acc_post = acc_post)
  
  return(res)
}