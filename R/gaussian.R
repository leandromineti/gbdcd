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
gbdcd <- function(y, neigh, c = 0.35, n_iterations = 1000000, burn_in = 500000,
                         mu0 = 0, sigma_0 = sqrt(2)) {

  # Normalization of the target variable
  mean_y <- mean(y)
  sd_y <- sd(y)
  y <- (y - mean_y) / sd_y

  n <- length(y)
  prob_dist <- (1 - c) ^ (1:n)
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

      new.partitions <- rcpp_partition(neigh, new_centers)
      new.means <- aggregate(y ~ new.partitions, FUN = mean)
      new.media <- subset(new.means, new.partitions == new_center)$y
      n.prop <- sum(new.partitions == new_center)

      mean.proposed <- (sigma2 / (n.prop * sigma2_0 + sigma2)) * mu0 +
        (n.prop * sigma2_0 / (n.prop * sigma2_0 + sigma2)) * new.media

      sigma.proposed <- sqrt(1 / ((1 / sigma2_0) + (n.prop / sigma2)))

      vec_means[new_center] <- rnorm(1, mean = mean.proposed, sd = sigma.proposed)

      phi <- dnorm(vec_means[new_center], mean = mean.proposed, sd = sigma.proposed)

      probMkplus1 <- dnorm(vec_means[new_center], mean = mu0, sd = sigma_0)

      # Likelihood ratio
      LLk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      LLkplus1 <- sum(dnorm(y, mean = vec_means[new.partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(LLkplus1 - LLk)
      k <- length(centers)

      A <- ratio * (1 - c) * 1 * (probMkplus1 / phi)

      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new.partitions
        v_accept[i] <- 1
      }
    }

    ## Death step -------------------------------------------------------------
    if (step == "Death") {
      center <- sample(centers, 1)
      means <- aggregate(y ~ partitions, FUN = mean)
      media <- subset(means, partitions == center)$y
      n.prop <- sum(partitions == center)

      mean.proposed <- (sigma2 / (n.prop * sigma2_0 + sigma2)) * mu0 +
        (n.prop * sigma2_0 / (n.prop * sigma2_0 + sigma2)) * media

      sigma.proposed <- sqrt(1 / ((1 / sigma2_0) + (n.prop / sigma2)))

      phi <- dnorm(vec_means[center], mean = mean.proposed, sd = sigma.proposed)

      probMkminus1 <- dnorm(vec_means[center], mean = mu0, sd = sigma_0)

      new_centers <- centers[-which(centers == center)]
      new.partitions <- rcpp_partition(neigh, new_centers)
      new.means <- aggregate(y ~ new.partitions, FUN = mean)

      # Likelihood ratio
      LLk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      LLkminus1 <- sum(dnorm(y, mean = vec_means[new.partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(LLkminus1 - LLk)
      k <- length(centers)

      A <- ratio * (1 / (1 - c)) * 1 * (phi / probMkminus1)

      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new.partitions
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

      new.partitions <- rcpp_partition(neigh, new_centers)

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
      LLk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      LLkshift <- sum(dnorm(y, mean = vec_means[new.partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(LLkshift - LLk)

      A <- ratio * (n_Gk / new.n_Gk) * (m_Gj / new.m_Gj)
      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new.partitions
        v_accept[i] <- 1
      }
    }

    ## Switch step ------------------------------------------------------------
    if (step == "Switch") {

      # Choose two centers and switch
      switch <- sample(1:length(centers), size = 2)
      new_centers <- replace(centers, list = switch, centers[rev(switch)])

      new.partitions <- rcpp_partition(neigh, new_centers)

      # Likelihood ratio
      LLk <- sum(dnorm(y, mean = vec_means[partitions], sd = sqrt(sigma2), log = TRUE))
      LLkswitch <- sum(dnorm(y, mean = vec_means[new.partitions], sd = sqrt(sigma2), log = TRUE))
      ratio <- exp(LLkswitch - LLk)

      A <- ratio
      alpha <- min(1, A)

      if (runif(1) < alpha) {
        centers <- new_centers
        partitions <- new.partitions
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
    if (i > burn_in) freq_matrix <- freq_matrix + RcppFreqMatrix(partitions)
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
