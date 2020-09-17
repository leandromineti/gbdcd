#' Gaussian BDCD
#'
#' Implementation of the Bayesian Detection of Clusters and Discontinuities
#'
#' @author leandromineti@gmail.com
#'
#' @param y a vector specifying the target variable for each of the nodes.
#' @param neigh a matrix defining the graph neighbors.
#' @param c parameter indicating the \emph{a priori} number of clusters.
#' @param n_iterations number of iterations.
#' @param burn_in number of discarded iterations.
#' @param mu0 \emph{a priori} mean.
#' @param sigma_0 \emph{a priori} standard deviation.
#'
#' @return a \code{list} of seven objects:
#' \itemize{
#'   \item mean.info: \emph{a posteriori} means and credible interval.
#'   \item cluster.info: hierarchical clustering return object.
#'   \item matConnections: frequency matrix indicating how many times each
#'   pair of nodes was in the same cluster.
#'   \item k.MCMC: a vector indicating the number of clusters in each iteration.
#'   \item mean_y: target variable mean.
#'   \item sd_y: target variable standard deviation.
#'   \item vec.centers: center configuration for each iteration.
#' }
#'
#' @export
#'
#' @examples
#' library(gbdcd)
#'
#' data("aneeldata", package = "gbdcd")
#' data("aneelshape", package = "gbdcd")
#'
#' target_variable <- aneelshape$z_Precipitation
#' neighbors <- aneeldata$connections
#'
#' out <- gbdcd(
#'   y = target_variable,
#'   neigh = neighbors,
#'   c = 0.35,
#'   n_iterations = 100000,
#'   burn_in = 50000,
#'   mu0 = 0,
#'   sigma_0 = sqrt(2)
#' )
#' @family gbdcd
gbdcd <- function(y,
                  neigh,
                  c = 0.35,
                  n_iterations = 1000000,
                  burn_in = 500000,
                  mu0 = 0,
                  sigma_0 = sqrt(2)) {

  # Normalization of the target variable
  mean_y <- mean(y)
  sd_y <- sd(y)
  y <- (y - mean_y) / sd_y

  n <- length(y)
  prob_dist <- (1 - c)^(1:n)
  prob_dist <- prob_dist / sum(prob_dist)
  mean_k <- round(sum((1:n) * prob_dist))

  sigma2_0 <- sigma_0^2 # A priori for the means in each group

  # Start-up variables
  k_vector <- rep(NA, n_iterations)
  v_sigma2 <- rep(NA, n_iterations)
  v_steps <- rep(NA, n_iterations)
  v_centers <- rep(NA, n_iterations)
  v_accept <- rep(0, n_iterations)
  mat_y_hat <- matrix(NA, n, n_iterations)
  freq_matrix <- matrix(0, n, n)
  vec_means <- rep(NA, n)
  centers <- sample.int(n, size = mean_k, replace = FALSE)
  partitions <- rcpp_partition(neigh, centers)

  # Initial cluster configuration
  means <- aggregate(y ~ partitions, FUN = mean)
  vec_means[means$partitions] <- means$y
  sigma2 <- sum((y - vec_means[partitions])^2) / (n - length(centers))

  # Initialize progress bar
  pb <- txtProgressBar(
    min = 0, max = n_iterations,
    char = "=", title = "progress bar"
  )

  ## Beginning of the chain ---------------------------------------------------
  for (i in 1:n_iterations) {
    # Chooses the step considering the present state of the clusters.
    k <- length(centers)

    if (k == n) {
      step <- sample(c("Death", "Update", "Switch"), size = 1, prob = c(0.8, 0.1, 0.1))
    } else {
      if (k == 1) {
        step <- sample(c("Birth", "Update", "Shift"), size = 1, prob = c(0.8, 0.1, 0.1))
      } else {
        step <- sample(c("Birth", "Death", "Update", "Shift", "Switch"),
          size = 1,
          prob = c(0.35, 0.35, 0.1, 0.1, 0.1)
        )
      }
    }

    ## Birth step -------------------------------------------------------------
    if (step == "Birth") {
      # Select a potential new center
      new_center <- sample((1:n)[!((1:n) %in% centers)], size = 1)
      # Insert the selected center into the vector of centers
      new_centers <- append(centers, new_center, after = sample(0:length(centers), size = 1))

      new_partitions <- rcpp_partition(neigh, new_centers)
      new_means <- aggregate(y ~ new_partitions, FUN = mean)
      new_mean <- subset(new_means, new_partitions == new_center)$y
      n_prop <- sum(new_partitions == new_center)

      mean_proposed <- (sigma2 / (n_prop * sigma2_0 + sigma2)) * mu0 +
        (n_prop * sigma2_0 / (n_prop * sigma2_0 + sigma2)) * new_mean

      sigma_proposed <- sqrt(1 / ((1 / sigma2_0) + (n_prop / sigma2)))

      vec_means[new_center] <- rnorm(1, mean = mean_proposed, sd = sigma_proposed)

      phi <- dnorm(vec_means[new_center], mean = mean_proposed, sd = sigma_proposed)

      probMkplus1 <- dnorm(vec_means[new_center], mean = mu0, sd = sigma_0)

      # Likelihood ratio
      llk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      llkplus1 <- sum(dnorm(y, mean = vec_means[new_partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(llkplus1 - llk)
      k <- length(centers)

      A <- ratio * (1 - c) * 1 * (probMkplus1 / phi)

      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new_partitions
        v_accept[i] <- 1
      }
    }

    ## Death step -------------------------------------------------------------
    if (step == "Death") {
      center <- sample(centers, 1)
      means <- aggregate(y ~ partitions, FUN = mean)
      media <- subset(means, partitions == center)$y
      n_prop <- sum(partitions == center)

      mean_proposed <- (sigma2 / (n_prop * sigma2_0 + sigma2)) * mu0 +
        (n_prop * sigma2_0 / (n_prop * sigma2_0 + sigma2)) * media

      sigma_proposed <- sqrt(1 / ((1 / sigma2_0) + (n_prop / sigma2)))

      phi <- dnorm(vec_means[center], mean = mean_proposed, sd = sigma_proposed)

      probMkminus1 <- dnorm(vec_means[center], mean = mu0, sd = sigma_0)

      new_centers <- centers[-which(centers == center)]
      new_partitions <- rcpp_partition(neigh, new_centers)
      new_means <- aggregate(y ~ new_partitions, FUN = mean)

      # Likelihood ratio
      llk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      llkminus1 <- sum(dnorm(y, mean = vec_means[new_partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(llkminus1 - llk)
      k <- length(centers)

      A <- ratio * (1 / (1 - c)) * 1 * (phi / probMkminus1)

      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new_partitions
        v_accept[i] <- 1
      }
    }

    ## Update step ------------------------------------------------------------
    if (step == "Update") {
      # Update the vector of means
      means <- aggregate(y ~ partitions, FUN = mean)
      count <- as.data.frame(table(partitions))
      means <- merge(means, count, by = "partitions")
      weigths <- sigma2 / (means$Freq * sigma2_0 + sigma2)

      # Update the vector os means: priori with local likelihood
      vec_means[means$partitions] <- rnorm(dim(means)[1],
        mean = weigths * mu0 + (1 - weigths) * means$y,
        sd = sqrt(1 / ((1 / sigma2_0) + (means$Freq / sigma2)))
      )

      # Update the variance parameter according to the priori and the likelihood
      # assuming known means
      k <- length(centers)

      if (k < n) {
        sum.2 <- sum((y - vec_means[partitions])^2)
        # Maximum likelihood
        S2 <- sum.2 / (n - k)
        S2.new <- geoR::rinvchisq(n = 1, df = n - k, scale = S2)

        if (S2.new < 100) {
          sigma2 <- S2.new
          v_accept[i] <- 1
        }
      }
      v_accept[i] <- 1
    }

    ## Shift step -------------------------------------------------------------
    if (step == "Shift") {
      centers_shift_ok <- centers
      new_centers <- centers

      for (cntr in centers) {
        # If all neighbors of a center are centers, remove the center from the vector
        if (all(neigh[which(neigh[, 1] == cntr), 2] %in% centers)) {
          centers_shift_ok <- centers_shift_ok[centers_shift_ok != cntr]
        }
      }

      # Given the possible centers, choose one
      if (length(centers_shift_ok) == 1) {
        center <- centers_shift_ok
      } else {
        center <- sample(centers_shift_ok, size = 1)
      }
      # Neighbors of the choosen center
      center_neighbors <- neigh[which(neigh[, 1] == center), 2]
      # Free neighbors of the choosen center
      free_neighbors <- center_neighbors[!(center_neighbors %in% centers)]

      if (length(free_neighbors) == 1) {
        shift <- free_neighbors
      } else {
        shift <- sample(free_neighbors, size = 1)
      }

      new_centers[which(centers == center)] <- shift

      n_Gk <- length(centers_shift_ok)
      m_Gj <- length(free_neighbors)

      new_partitions <- rcpp_partition(neigh, new_centers)

      vec_means[shift] <- vec_means[center]

      new_centers_shift_ok <- new_centers
      # Repeats the previous procedure for the new vector of centers
      for (nc in new_centers) {
        if (all(neigh[which(neigh[, 1] == nc), 2] %in% new_centers)) {
          new_centers_shift_ok <- new_centers_shift_ok[new_centers_shift_ok != nc]
        }
      }

      # Neighbors of the choosen center
      neigh_new_center <- neigh[which(neigh[, 1] == shift), 2]
      # Free neighbors of the choosen center
      new.free_neighbors <- neigh_new_center[!(neigh_new_center %in% new_centers)]

      new.n_Gk <- length(new_centers_shift_ok)
      new.m_Gj <- length(new.free_neighbors)

      # Likelihood ratio
      llk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      llkshift <- sum(dnorm(y, mean = vec_means[new_partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(llkshift - llk)

      A <- ratio * (n_Gk / new.n_Gk) * (m_Gj / new.m_Gj)
      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new_partitions
        v_accept[i] <- 1
      }
    }

    ## Switch step ------------------------------------------------------------
    if (step == "Switch") {

      # Choose two centers and switch
      switch <- sample(1:length(centers), size = 2)
      new_centers <- replace(centers, list = switch, centers[rev(switch)])

      new_partitions <- rcpp_partition(neigh, new_centers)

      # Likelihood ratio
      llk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      llk_switch <- sum(dnorm(y, mean = vec_means[new_partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(llk_switch - llk)

      A <- ratio
      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new_partitions
        v_accept[i] <- 1
      }
    }

    # Update output variables
    k_vector[i] <- length(centers)
    v_sigma2[i] <- sigma2
    v_steps[i] <- step
    mat_y_hat[, i] <- vec_means[partitions]
    v_centers[i] <- paste(centers, collapse = ";")

    # Set progress bar
    setTxtProgressBar(pb, i, "running Reversible Jump...")

    # Frequency matrix data
    if (i > burn_in) freq_matrix <- freq_matrix + rcpp_freqmatrix(partitions)
  }

  Sys.sleep(1)
  close(pb) # Ends progress bar

  ## Results processing -------------------------------------------------------

  # Remove burn-in
  seq.burn <- -(1:burn_in)
  k_vector <- k_vector[seq.burn]
  v_sigma2 <- v_sigma2[seq.burn]
  v_steps <- v_steps[seq.burn]
  mat_y_hat <- mat_y_hat[, seq.burn]
  v_accept <- v_accept[seq.burn]
  v_centers <- v_centers[seq.burn]

  # A posteriori means and confidence interval
  Yhat <- apply(mat_y_hat, MARGIN = 1, FUN = median)
  lwr <- apply(mat_y_hat, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.05))
  upr <- apply(mat_y_hat, MARGIN = 1, FUN = function(x) quantile(x, probs = 0.95))

  # Processing frequency matrix results
  matConnections <- freq_matrix
  matConnections[lower.tri(matConnections)] <- t(matConnections)[lower.tri(matConnections)]

  maximum <- max(as.vector(matConnections), na.rm = TRUE)
  matConnections_mod <- maximum - matConnections
  diag(matConnections_mod) <- maximum + 1

  # Partitioning using hierarchical clustering
  clusters <- hclust(as.dist(matConnections_mod), c("single", "complete")[1])

  output <- list(
    mean.info = cbind(lwr, Yhat, upr),
    cluster.info = clusters,
    matConnections = matConnections,
    k.MCMC = k_vector,
    mean_y = mean_y,
    sd_y = sd_y,
    vec.centers = v_centers
  )

  return(output)
}

#' Gaussian BDCD regression
#'
#' Implementation of the Bayesian Detection of Clusters and Discontinuities
#'
#' @author leandromineti@gmail.com
#'
#' @param y a vector specifying the target variable for each of the nodes.
#' @param x a matrix with the independent variables.
#' @param viz a matrix defining the graph neighbors.
#' @param n_iterations number of iterations.
#' @param burn_in number of discarded iterations.
#' @param c parameter indicating the \emph{a priori} number of clusters.
#' @param prior_coeffs_mu \emph{a priori} mean for regresion coefficients.
#' @param prior_sigma \emph{a priori} standard deviation.
#' @param lambda regularization parameter.
#' @param plot plot the results. Default = FALSE.
#'
#' @return a \code{list} of seven objects:
#' \itemize{
#'   \item mean.info: \emph{a posteriori} means and credible interval.
#'   \item cluster.info: hierarchical clustering return object.
#'   \item matConnections: frequency matrix indicating how many times each
#'   pair of nodes was in the same cluster.
#'   \item k.MCMC: a vector indicating the number of clusters in each iteration.
#'   \item mean_y: target variable mean.
#'   \item sd_y: target variable standard deviation.
#'   \item vec.centers: center configuration for each iteration.
#' }
#'
#' @export
#'
#' @family gbdcd
gbdcd_regression <- function(y,
                             x,
                             viz,
                             n_iterations = 100000,
                             burn_in = 50000,
                             c = 0.35,
                             prior_coeffs_mu = rep(0, 2),
                             prior_sigma = sqrt(2),
                             lambda = 0.001,
                             plot = F) {

  # Define important simulation parameters
  n_regions <- length(y)
  sigma2 <- prior_sigma^2
  a_0 <- 2.1 # Inverse Gamma shape parameter
  b_0 <- 1.1 # Inverse Gamma scale parameter

  n_coeffs <- dim(x)[2] # Number of estimated coefficients
  n_coeffs_names <- paste0("coeff_", 0:(n_coeffs - 1)) # Coefficient names

  # Covariancia a priori dos betas (falta multiplicar pelo sigma2)
  inv_lambda_identity_prior <- diag(rep(1 / lambda, 2))

  # Bookkeeping variables
  vec_k <- rep(NA, n_iterations)
  vec_sigma2 <- rep(NA, n_iterations)
  vec_steps <- rep(NA, n_iterations)
  vec_centers <- rep(NA, n_iterations)
  vec_accept <- rep(0, n_iterations)
  mat_freq <- matrix(0, n_regions, n_regions)
  mat_coeffs <- matrix(NA, n_regions, ncol = n_coeffs)
  mat_coeffs_hat <- array(NA, dim = c(n_regions, n_coeffs, n_iterations), dimnames = list(NULL, n_coeffs_names, NULL))
  cluster_centers <- sample.int(n_regions, size = 2, replace = FALSE)
  cluster_partitions <- rcpp_partition(viz, cluster_centers)

  # Inicializa com os mesmos betas para todos os clusters
  mat_coeffs[cluster_centers, ] <- 0

  # Start execution bar
  pb <- txtProgressBar(min = 0, max = n_iterations, char = "=", title = "progress bar")

  # Start the Markov Chain Monte Carlo simulation
  for (i in 1:n_iterations)
  {
    k <- length(cluster_centers) # Number of cluster centers on the current state

    # Select step type
    if (k == 1) {
      step_choice <- sample(c("birth", "death", "update"), size = 1, prob = c(0.8, 0, 0.2))
    } else if (k == n_regions) {
      step_choice <- sample(c("birth", "death", "update"), size = 1, prob = c(0, 0.8, 0.2))
    } else {
      step_choice <- sample(c("birth", "death", "update"), size = 1, prob = c(0.4, 0.4, 0.2))
    }

    # Define birth step
    if (step_choice == "birth") {

      # Select potential new cluster center and update cluster configuration
      new_center <- sample((1:n_regions)[!((1:n_regions) %in% cluster_centers)], size = 1)
      new_cluster_centers <- append(cluster_centers, new_center,
        after = sample(0:length(cluster_centers), size = 1)
      )
      new_cluster_partitions <- rcpp_partition(viz, new_cluster_centers)
      y_cluster <- subset(y, new_cluster_partitions == new_center)
      x_cluster <- subset(x, new_cluster_partitions == new_center)

      # Propose regression coefficients for the new cluster
      auxiliar_proposed <- solve(crossprod(x_cluster, x_cluster) + lambda * diag(2))
      coeffs_mean_proposed <- auxiliar_proposed %*% crossprod(x_cluster, y_cluster)
      coeffs_cov_proposed <- sigma2 * auxiliar_proposed

      mat_coeffs[new_center, ] <- mvtnorm::rmvnorm(1,
        mean = coeffs_mean_proposed,
        sigma = coeffs_cov_proposed
      )

      phi <- mvtnorm::dmvnorm(mat_coeffs[new_center, ],
        mean = coeffs_mean_proposed,
        sigma = coeffs_cov_proposed,
        log = TRUE
      )

      phi_k_plus_1 <- mvtnorm::dmvnorm(mat_coeffs[new_center, ],
        mean = prior_coeffs_mu,
        sigma = sigma2 * inv_lambda_identity_prior,
        log = TRUE
      )

      # Define likelihood ratio and step acceptance probability
      llk <- sum(dnorm(y,
        mean = rowSums(mat_coeffs[cluster_partitions, ] * x),
        sd = sqrt(sigma2), log = TRUE
      ))

      llk_plus_1 <- sum(dnorm(y,
        mean = rowSums(mat_coeffs[new_cluster_partitions, ] * x),
        sd = sqrt(sigma2), log = TRUE
      ))

      llk_ratio <- exp(llk_plus_1 - llk)
      accept_prob <- llk_ratio * (1 - c) * 1 * exp(phi_k_plus_1 - phi)

      if (is.na(accept_prob)) {
        accept_prob <- 0
      } # Handle errors on small probabilities

      alpha <- min(1, accept_prob)

      if (runif(1) < alpha) {
        cluster_centers <- new_cluster_centers
        cluster_partitions <- new_cluster_partitions
        vec_accept[i] <- 1
      }
    }

    # Define death step
    if (step_choice == "death") {

      # Select potential cluster center to remove and update cluster configuration
      center <- sample(cluster_centers, 1)
      new_cluster_centers <- cluster_centers[-which(cluster_centers == center)]
      new_cluster_partitions <- rcpp_partition(viz, new_cluster_centers)
      y_cluster <- subset(y, cluster_partitions == center)
      x_cluster <- subset(x, cluster_partitions == center)

      # Calculate coefficient distribution parameters for the proposed configuration
      auxiliar_proposed <- solve(crossprod(x_cluster, x_cluster) + lambda * diag(2))
      coeffs_mean_proposed <- auxiliar_proposed %*% crossprod(x_cluster, y_cluster)
      coeffs_cov_proposed <- sigma2 * auxiliar_proposed

      phi <- mvtnorm::dmvnorm(mat_coeffs[center, ],
        mean = coeffs_mean_proposed,
        sigma = coeffs_cov_proposed,
        log = TRUE
      )

      phi_k_minus_1 <- mvtnorm::dmvnorm(mat_coeffs[center, ],
        mean = prior_coeffs_mu,
        sigma = sigma2 * inv_lambda_identity_prior,
        log = TRUE
      )

      # Define likelihood ratio and step acceptance probability
      llk <- sum(dnorm(y,
        mean = rowSums(mat_coeffs[cluster_partitions, ] * x),
        sd = sqrt(sigma2), log = TRUE
      ))

      llk_minus_1 <- sum(dnorm(y,
        mean = rowSums(mat_coeffs[new_cluster_partitions, ] * x),
        sd = sqrt(sigma2), log = TRUE
      ))

      llk_ratio <- exp(llk_minus_1 - llk)
      accept_prob <- llk_ratio * (1 / (1 - c)) * 1 * exp(phi - phi_k_minus_1)

      if (is.na(accept_prob)) {
        accept_prob <- 0
      } # Handle errors on small probabilities

      alpha <- min(1, accept_prob)

      if (runif(1) < alpha) {
        cluster_centers <- new_cluster_centers
        cluster_partitions <- new_cluster_partitions
        vec_accept[i] <- 1
      }
    }

    # Define update step
    if (step_choice == "update") {

      # Update regression coefficients
      for (center_update in cluster_centers) {

        # Select the data points of an specific cluster
        y_cluster <- subset(y, cluster_partitions == center_update)
        x_cluster <- subset(x, cluster_partitions == center_update)

        # Update regression coefficients on an specific cluster
        auxiliar_proposed <- solve(crossprod(x_cluster, x_cluster) + lambda * diag(2))
        coeffs_mean_proposed <- auxiliar_proposed %*% crossprod(x_cluster, y_cluster)
        coeffs_cov_proposed <- sigma2 * auxiliar_proposed

        mat_coeffs[center_update, ] <- mvtnorm::rmvnorm(1,
          mean = coeffs_mean_proposed,
          sigma = coeffs_cov_proposed
        )
      }

      # Update variance parameter
      a_n <- a_0 + n_regions / 2
      b_n <- b_0 + sum((y - rowSums(mat_coeffs[cluster_partitions, ] * x))^2) / 2
      sigma2 <- 1 / rgamma(1, shape = a_n, rate = b_n)

      vec_accept[i] <- 1
    }

    # Record chain state
    vec_k[i] <- length(cluster_centers)
    vec_sigma2[i] <- sigma2
    vec_steps[i] <- step_choice
    mat_coeffs_hat[, , i] <- mat_coeffs[cluster_partitions, ]
    vec_centers[i] <- paste(cluster_centers, collapse = ";")

    # Set progress bar
    setTxtProgressBar(pb, i)

    # Fill 'connection' frequency matrix
    if (i > burn_in) mat_freq <- mat_freq + rcpp_freqmatrix(cluster_partitions)
  }

  # End of MCMC iterations
  Sys.sleep(1)
  close(pb)

  # Chain results processing - remove burn-in
  seq.burn <- -(1:burn_in)
  vec_k <- vec_k[seq.burn]
  vec_sigma2 <- vec_sigma2[seq.burn]
  vec_steps <- vec_steps[seq.burn]
  mat_coeffs_hat <- mat_coeffs_hat[, , seq.burn]
  vec_accept <- vec_accept[seq.burn]
  vec_centers <- vec_centers[seq.burn]

  # Processing 'connection' frequency matrix
  mat_connections <- mat_freq
  mat_connections[lower.tri(mat_connections)] <- t(mat_connections)[lower.tri(mat_connections)]
  maximum <- max(as.vector(mat_connections), na.rm = TRUE)
  mat_connections <- maximum - mat_connections
  diag(mat_connections) <- maximum + 1
  clusters <- hclust(as.dist(mat_connections), c("single", "complete")[1])

  # Plot results
  if (plot == TRUE) {
    plot(1:length(vec_k), vec_k, type = "l", xlab = "step", main = "k-MCMC")

    barplot(table(vec_k) / sum(table(vec_k)),
      col = "light blue",
      main = "k-Posteriori"
    )

    cat("\n Estimates for k: \n")
    print(summary(vec_k))
    print(quantile(vec_k, probs = c(0.05, 0.95)))

    cat("Step acceptance frequency: \n \n")
    print(aggregate(vec_accept ~ vec_steps, FUN = mean))

    cat("\n Estimates for sigma2: \n")
    print(summary(vec_sigma2))
    print(quantile(vec_sigma2, probs = c(0.05, 0.95)))

    plot(clusters,
      ylab = "proximity",
      main = "Proximity Dendogram",
      xlab = "", cex = 0.7
    )
  }

  output <- list(
    coeff.info = mat_coeffs_hat,
    cluster.info = clusters,
    matConexoes = mat_connections,
    k.MCMC = vec_k,
    vec.centros = vec_centers,
    vec.sigma2 = vec_sigma2
  )

  return(output)
}

#' GBDCD groups
#'
#' Implement Adapted Ng-Jordan-Weiss clustering on gbdcd results
#'
#' @author leandromineti@gmail.com
#'
#' @param out output from gbdcd.
#' @param k number of clusters
#'
#' @return a \code{vector} indicating the cluster for each element.
#'
#' @export
#'
#' @family gbdcd
gbdcdGroups <- function(out, k = 2) {
  S <- out$matConnections
  S <- S / length(out$vec.centers)

  # Appendix B: Adapted Ng-Jordan-Weiss algorithm

  diag_aux <- diag(1 / sqrt(rowSums(S)))
  L <- diag_aux %*% S %*% diag_aux
  # Getting eigenvalues and eigenvectors
  eigen_output <- eigen(L)
  eigen_results <- eigen_output$vectors[, order(eigen_output$values, decreasing = T)[1:k]]

  fit <- kmeans(eigen_results, k)

  return(fit$cluster)
}

# Clean up when the package is unloaded
.onUnload <- function(libpath) {
  library.dynam.unload("gbdcd", libpath)
}

# Some needed roxygen2 tags
#' @useDynLib gbdcd
#' @importFrom Rcpp sourceCpp
NULL
# > NULL
