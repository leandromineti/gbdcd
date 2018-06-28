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
#' @param sigma0 \emph{a priori} standard deviation.
#'
#' @return a \code{list} of seven objects:
#' \itemize{
#'   \item mean.info: \emph{a posteriori} means and credible interval.
#'   \item cluster.info: hierarchical clustering return object.
#'   \item matConnections: frequency matrix indicating how many times each pair of
#'   nodes was in the same cluster.
#'   \item k.MCMC: a vector indicating the number of clusters in each iteration.
#'   \item mean.y: target variable mean.
#'   \item sd.y: target variable standard deviation.
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
#' out <- gaussianBDCD(y = target_variable, 
#'                    neigh = neighbors, 
#'                    c = 0.35, 
#'                    n_iterations = 100000, 
#'                    burn_in = 50000, 
#'                    mu0 = 0, 
#'                    sigma0 = sqrt(2))
#'
#' @family gbdcd
gaussianBDCD <- function(y, neigh, c = 0.35, n_iterations = 1000000, burn_in = 500000,
                         mu0=0, sigma0=sqrt(2)) {
  
  # Normalization of the target variable
  mean.y <- mean(y)
  sd.y <- sd(y)
  y <- (y - mean.y)/sd.y
  
  N <- length(y)
  Probs <- (1-c)^(1:N)
  Probs <- Probs/sum(Probs)
  mean.k <- round(sum((1:N)*Probs))  # A priori mean for the number of clusters.
  
  sig2ma0 <- sigma0^2 # A priori for the means in each group
  
  # Start-up variables
  k_vector <- rep(NA, n_iterations)
  v_sigma2 <- rep(NA, n_iterations)
  v_steps <- rep(NA, n_iterations)
  v_centers <- rep(NA, n_iterations)
  v_accept <- rep(0, n_iterations)
  mat.Yhat <- matrix(NA, N, n_iterations)
  freq_matrix <- matrix(0, N, N)
  vec.means <- rep(NA, N)
  centers <- sample.int(N, size = mean.k, replace = FALSE) 
  partitions <- RcppPartition(neigh, centers)
  
  # Initial cluster configuration
  means <- aggregate(y ~ partitions, FUN=mean)
  vec.means[means$partitions] <- means$y
  sigma2 <- sum( (y - vec.means[partitions])^2 )/(N-length(centers))
  
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iterations, char= "=", title= "progress bar")
  
  ## Beginning of the chain ---------------------------------------------------
  for(i in 1:n_iterations) {
    # Chooses the step considering the present state of the clusters.
    k <- length(centers)
    
    if(k == N) {
      step <- sample(c("Death", "Update", "Switch"), size=1, prob = c(0.8, 0.1, 0.1))
    } else { 
      if (k == 1) {
        step <- sample(c("Birth", "Update", "Shift"), size=1, prob = c(0.8, 0.1, 0.1))
      } else {
        step <- sample(c("Birth", "Death", "Update", "Shift", "Switch"), size=1, 
                        prob = c(0.35, 0.35, 0.1, 0.1, 0.1))
      }
    }
    
    ## Birth step -------------------------------------------------------------
    if(step == "Birth") {
      # Select a potential new center
      new.center  <- sample((1:N)[!((1:N) %in% centers)], size=1)
      # Insert the selected center into the vector of centers
      new.centers <- append(centers, new.center, after = sample(0:length(centers), size=1))
      
      new.partitions <- RcppPartition(neigh, new.centers)
      new.means <- aggregate(y ~ new.partitions, FUN=mean)
      new.media <- subset(new.means,  new.partitions == new.center)$y
      n.prop <- sum(new.partitions == new.center)
      
      mean.proposed  <- (sigma2/(n.prop*sig2ma0 + sigma2))*mu0 + 
        (n.prop*sig2ma0/(n.prop*sig2ma0 + sigma2))*new.media
      
      sigma.proposed  <- sqrt(1/( (1/sig2ma0) +  (n.prop/sigma2)))
      
      vec.means[new.center] <- rnorm(1, mean = mean.proposed, sd = sigma.proposed)
      
      phi <- dnorm(vec.means[new.center], mean = mean.proposed, sd = sigma.proposed)
      
      probMkplus1 <- dnorm(vec.means[new.center], mean = mu0, sd = sigma0)
      
      # Likelihood ratio
      LLk <- sum(dnorm(y, mean = vec.means[partitions], sd=sqrt(sigma2), log=TRUE))
      LLkplus1 <- sum(dnorm(y, mean = vec.means[new.partitions], sd=sqrt(sigma2), log=TRUE))
      ratio <- exp(LLkplus1 - LLk)
      k <- length(centers)
      
      A <- ratio * (1-c) * 1 * (probMkplus1/phi)
      
      alpha <- min(1, A)
      
      if(runif(1) < alpha) { 
        centers <- new.centers
        partitions <- new.partitions
        v_accept[i] <- 1
      }
    }

    ## Death step -------------------------------------------------------------
    if(step == "Death") { 
      center <- sample(centers, 1)
      means <- aggregate(y ~ partitions, FUN=mean)
      media <- subset(means,  partitions == center)$y
      n.prop <- sum(partitions == center)
      
      mean.proposed <- (sigma2/(n.prop*sig2ma0 + sigma2))*mu0 + 
        (n.prop*sig2ma0/(n.prop*sig2ma0 + sigma2))*media
      
      sigma.proposed  <- sqrt(1/( (1/sig2ma0) +  (n.prop/sigma2)))
      
      phi <- dnorm(vec.means[center], mean = mean.proposed, sd = sigma.proposed) 
      
      probMkminus1 <- dnorm(vec.means[center], mean = mu0, sd = sigma0) 
      
      new.centers <- centers[-which(centers==center)] 
      new.partitions <- RcppPartition(neigh, new.centers)
      new.means <- aggregate(y ~ new.partitions, FUN=mean)
      
      # Likelihood ratio
      LLk <- sum(dnorm(y, mean = vec.means[partitions],     sd=sqrt(sigma2), log=TRUE))
      LLkminus1 <- sum(dnorm(y, mean = vec.means[new.partitions], sd=sqrt(sigma2), log=TRUE))
      ratio <- exp(LLkminus1 - LLk)
      k <- length(centers)

      A <- ratio * (1/(1-c)) * 1 * (phi/probMkminus1)
      
      alpha <- min(1, A)
      
      if(runif(1) < alpha) { 
        centers <- new.centers
        partitions <- new.partitions
        v_accept[i] <- 1
      } 
    }
    
    ## Update step ------------------------------------------------------------
    if(step == "Update") {
      # Update the vector of means
      means <- aggregate(y ~ partitions, FUN=mean)
      count <- as.data.frame(table(partitions))
      means <- merge(means, count, by="partitions") 
      weigths <- sigma2/(means$Freq*sig2ma0 + sigma2)
      
      # Update the vector os means: priori with local likelihood
      vec.means[means$partitions] <- rnorm(dim(means)[1],
                                            mean = weigths*mu0 + (1-weigths)*means$y,
                                            sd   = sqrt( 1/((1/sig2ma0) + (means$Freq/sigma2)) )  )
      
      # Update the variance parameter according to the priori and the likelihood
      # assuming known means
      k <- length(centers)
      
      if(k < N) {
        sum.2  <- sum( (y - vec.means[partitions])^2 ) 
        # Maximum likelihood
        S2 <- sum.2/(N-k)
        S2.new <- geoR::rinvchisq(n=1, df=N-k, scale=S2) 
        
        if(S2.new < 100) {
          sigma2 <- S2.new
          v_accept[i] <- 1
        }
      }
      v_accept[i] <- 1
    }
    
    ## Shift step -------------------------------------------------------------
    if(step == "Shift") {
      
      centers_shift_ok <- centers
      new.centers <- centers
      
      for(cntr in centers) {
        # If all neighbors of a center are centers, remove the center from the vector
        if(all(neigh[which(neigh[,1]==cntr),2] %in% centers)) 
          centers_shift_ok <- centers_shift_ok[centers_shift_ok!=cntr]
      }
      
      # Given the possible centers, choose one
      if(length(centers_shift_ok) == 1) { 
        center <- centers_shift_ok
      } else   center <- sample(centers_shift_ok, size = 1)
      # Neighbors of the choosen center
      center_neighbors <- neigh[which(neigh[,1]==center),2]
      # Free neighbors of the choosen center
      free_neighbors  <- center_neighbors[!(center_neighbors %in% centers)]
      
      if(length(free_neighbors)==1) { 
        shift <- free_neighbors
      } else shift <- sample(free_neighbors, size = 1)
      
      new.centers[which(centers==center)] <- shift
      
      n_Gk <- length(centers_shift_ok)
      m_Gj <- length(free_neighbors)
      
      new.partitions  <- RcppPartition(neigh, new.centers)
      
      vec.means[shift] <- vec.means[center]
      
      new.centers_shift_ok <- new.centers
      # Repeats the previous procedure for the new vector of centers
      for(nc in new.centers) {
        if(all(neigh[which(neigh[,1]==nc),2] %in% new.centers)) 
          new.centers_shift_ok <- new.centers_shift_ok[new.centers_shift_ok!=nc]
      }
      
      # Neighbors of the choosen center
      neigh_new_center <- neigh[which(neigh[,1]==shift),2]
      # Free neighbors of the choosen center
      new.free_neighbors <- neigh_new_center[!(neigh_new_center %in% new.centers)]
      
      new.n_Gk <- length(new.centers_shift_ok)
      new.m_Gj <- length(new.free_neighbors)
      
      # Likelihood ratio
      LLk <- sum(dnorm(y, mean = vec.means[partitions],     sd=sqrt(sigma2), log=TRUE))
      LLkshift <- sum(dnorm(y, mean = vec.means[new.partitions], sd=sqrt(sigma2), log=TRUE))
      ratio <- exp(LLkshift - LLk)
      
      A = ratio * (n_Gk/new.n_Gk) * (m_Gj/new.m_Gj)
      alpha <- min(1, A)
      
      if(runif(1) < alpha) {
        centers <- new.centers
        partitions <- new.partitions
        v_accept[i] <- 1
      }
    }
    
    ## Switch step ------------------------------------------------------------
    if(step == "Switch") {
      
      # Choose two centers and switch
      switch <- sample(1:length(centers), size=2)
      new.centers <- replace(centers, list = switch, centers[rev(switch)])
      
      new.partitions <- RcppPartition(neigh, new.centers)
      
      # Likelihood ratio
      LLk <- sum(dnorm(y, mean = vec.means[partitions],     sd=sqrt(sigma2), log=TRUE))
      LLkswitch <- sum(dnorm(y, mean = vec.means[new.partitions], sd=sqrt(sigma2), log=TRUE))
      ratio <- exp(LLkswitch - LLk)
      
      A = ratio
      alpha <- min(1, A)
      
      if(runif(1) < alpha) {
        centers <- new.centers
        partitions <- new.partitions
        v_accept[i] <- 1
      }
    }
    
    # Update output variables
    k_vector[i] <- length(centers)
    v_sigma2[i] <- sigma2
    v_steps[i] <- step
    mat.Yhat[,i] <- vec.means[partitions]
    v_centers[i] <- paste(centers, collapse=";")
    
    # Set progress bar	
    setTxtProgressBar(pb, i, "running Reversible Jump...")
    
    # Frequency matrix data
    if(i > burn_in) freq_matrix <- freq_matrix + RcppFreqMatrix(partitions)
  }
  
  Sys.sleep(1)
  close(pb)  # Ends progress bar
  
  ## Results processing -------------------------------------------------------
  
  # Remove burn-in
  seq.burn <- -(1:burn_in)
  k_vector <- k_vector[seq.burn]
  v_sigma2 <- v_sigma2[seq.burn]
  v_steps <- v_steps[seq.burn]
  mat.Yhat <- mat.Yhat[,seq.burn]
  v_accept <- v_accept[seq.burn]
  v_centers <- v_centers[seq.burn]
  
  # A posteriori means and confidence interval
  Yhat <- apply(mat.Yhat, MARGIN=1, FUN=median)
  lwr <- apply(mat.Yhat, MARGIN=1, FUN=function(x) quantile(x, probs=0.05))
  upr <- apply(mat.Yhat, MARGIN=1, FUN=function(x) quantile(x, probs=0.95))
  
  # Processing frequency matrix results
  matConnections <- freq_matrix
  matConnections[lower.tri(matConnections)] <- t(matConnections)[lower.tri(matConnections)]
  
  maximum <- max(as.vector(matConnections), na.rm=TRUE)
  matConnections_mod <- maximum - matConnections
  diag(matConnections_mod) <- maximum + 1
  
  # Partitioning using hierarchical clustering 
  clusters <- hclust(as.dist(matConnections_mod), c("single","complete")[1])
  
  output <- list(mean.info    = cbind(lwr, Yhat, upr),
                cluster.info = clusters,
                matConnections  = matConnections,
                k.MCMC       = k_vector,
                mean.y       = mean.y,
                sd.y         = sd.y,
                vec.centers  = v_centers)
  
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
#' @return a \code{list} of seven objects:
#' \itemize{
#'   \item mydata: k eigenvectors.
#'   \item fit: result of k-means on eigenvectors. 
#'   \item groups: grouped elements.
#' }
#'
#' @export
#'
#' @family gbdcd
gbdcdGroups <- function(out, k=2) {
  S <- out$matConnections
  S <- S/length(out$vec.centers)
  
  ## APPENDIX B: Adapted Ng-Jordan-Weiss algorithm
  D <- diag(rowSums(S))
  Daux <- diag(1/sqrt(rowSums(S)))
  L <- Daux%*%S%*%Daux
  ## Getting eigenvalues and eigenvectors
  saida <- eigen(L)
  vetores <- saida$vectors[, order(saida$values,decreasing=T)[1:k]]
  
  output <- list(mydata = vetores, 
                 fit = kmeans(mydata, k),
                 groups = fit$cluster)
  
  return(output)
}
