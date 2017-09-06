#' Gaussian BDCD
#' 
#' Implementation of the Bayesian Detection of Clusters and Discontinuities
#'
#' @author leandromineti@gmail.com
#'
#' @param y a vector of indices.
#' @param viz a matrix of neighbors.
#' @param c parameter indicating the a priori number of clusters.
#' @param n_iterations number of iterations.
#' @param burn_in number of discarded iterations.
#' @param mu0 priori.
#' @param sigma0 priori.
#'
#' @return \code{list} 
#'
#' @export
#' 
#' @examples
#'
#' @family gbdcd
gaussianBDCD <- function(y, viz, c = 0.35, n_iterations = 1000000, burn_in = 500000,
                         mu0=0, sigma0=sqrt(2)) {
  
  # Normalization of the target variable
  mean.y <- mean(y)
  sd.y <- sd(y)
  y <- (y - mean.y)/sd.y
  
  N <- length(y)
  Probs <- (1-c)^(1:N)
  Probs <- Probs/sum(Probs)
  media.k <- round( sum((1:N)*Probs) )  # A priori mean for the number of clusters.
  
  # A priori for the means in each group
  sig2ma0   <- sigma0^2   
  
  # Start-up variables
  k_vector       <- rep(NA, n_iterations)
  v_sigma2       <- rep(NA, n_iterations)
  v_passos       <- rep(NA, n_iterations)
  v_centros      <- rep(NA, n_iterations)
  v_aceite       <- rep(0, n_iterations)
  mat.Yhat       <- matrix(NA, N, n_iterations)
  freq_matrix    <- matrix(0, N, N)
  vec.Medias     <- rep(NA, N)
  centros        <- sample.int(N, size = media.k, 
                               replace = FALSE) 
  particoes      <- RcppPartition(viz, centros)
  
  # Initial cluster configuration
  medias <- aggregate(y ~ particoes, FUN=mean)
  vec.Medias[medias$particoes] <- medias$y
  sigma2 <- sum( (y - vec.Medias[particoes])^2 )/(N-length(centros))
  
  # Initializar progress bar
  pb <- txtProgressBar(min = 0, max = n_iterations, char= "=", title= "progress bar")
  
  ## Beginning of the chain ---------------------------------------------------
  for(i in 1:n_iterations)
  {
    # Chooses the step considering the present state of the clusters.
    k <- length(centros)
    
    if(k == N){
      passo <- sample(c("Death", "Update", "Switch"), size=1, 
                      prob = c(0.8, 0.1, 0.1))
    } else { 
      if (k == 1){
        passo <- sample(c("Birth", "Update", "Shift"), size=1, 
                        prob = c(0.8, 0.1, 0.1))
      } else {
        passo <- sample(c("Birth", "Death", "Update", "Shift", "Switch"), size=1, 
                        prob = c(0.35, 0.35, 0.1, 0.1, 0.1))
      }
    }
    
    ## Birth step -------------------------------------------------------------
    if(passo == "Birth"){
      ## Seleciona um potencial novo centro
      new.centro  <- sample((1:N)[!((1:N) %in% centros)], size=1)
      ## Insere no vetor de centros
      new.centros <- append(centros, new.centro, after = sample(0:length(centros), size=1))
      
      new.particoes   <- RcppPartition(viz, new.centros)
      new.medias      <- aggregate(y ~ new.particoes, FUN=mean)
      new.media       <- subset(new.medias,  new.particoes == new.centro)$y
      n.prop          <- sum(new.particoes == new.centro)
      
      media.proposta  <- (sigma2/(n.prop*sig2ma0 + sigma2))*mu0 + 
        (n.prop*sig2ma0/(n.prop*sig2ma0 + sigma2))*new.media
      
      sigma.proposta  <- sqrt(1/( (1/sig2ma0) +  (n.prop/sigma2)))
      
      vec.Medias[new.centro] <- rnorm(1, mean = media.proposta, 
                                      sd = sigma.proposta)
      
      phi <- dnorm(vec.Medias[new.centro], mean = media.proposta, 
                   sd = sigma.proposta)
      
      probMkplus1 <- dnorm(vec.Medias[new.centro], mean = mu0, 
                           sd = sigma0)
      
      ## Cálculo da razão de Verossimilhança
      LLk      <- sum( dnorm(y, mean = vec.Medias[particoes],     sd=sqrt(sigma2), log=TRUE ) )
      LLkplus1 <- sum( dnorm(y, mean = vec.Medias[new.particoes], sd=sqrt(sigma2), log=TRUE ) )
      RAZAO    <- exp(LLkplus1 - LLk)
      k        <- length(centros)
      
      ## A     <- RAZÃO * (1-c) * (1/(N-k)) * 1 * (probMkplus1/phi) ## ORIGINAL
      A     <- RAZAO * (1-c) * 1 * (probMkplus1/phi)
      
      alpha <- min(1, A)
      
      if(runif(1) < alpha){ 
        centros   <- new.centros
        particoes <- new.particoes
        v_aceite[i] <- 1
      }
    }

    ## Death step -------------------------------------------------------------
    if(passo == "Death"){ 
      centro  <- sample(centros, 1)
      medias  <- aggregate(y ~ particoes, FUN=mean)
      media   <- subset(medias,  particoes == centro)$y
      n.prop  <- sum(particoes == centro)
      
      media.proposta  <- (sigma2/(n.prop*sig2ma0 + sigma2))*mu0 + 
        (n.prop*sig2ma0/(n.prop*sig2ma0 + sigma2))*media
      
      sigma.proposta  <- sqrt(1/( (1/sig2ma0) +  (n.prop/sigma2)))
      
      phi     <- dnorm(vec.Medias[centro], mean = media.proposta, 
                       sd = sigma.proposta) 
      
      probMkminus1 <- dnorm(vec.Medias[centro], mean = mu0, 
                            sd = sigma0) 
      
      new.centros   <- centros[-which(centros==centro)] 
      new.particoes <- RcppPartition(viz, new.centros)
      new.medias    <- aggregate(y ~ new.particoes, FUN=mean)
      
      ## Calculo da razão de Verossimilhança
      LLk       <- sum( dnorm(y, mean = vec.Medias[particoes],     sd=sqrt(sigma2), log=TRUE ) )
      LLkminus1 <- sum( dnorm(y, mean = vec.Medias[new.particoes], sd=sqrt(sigma2), log=TRUE ) )
      RAZAO     <- exp(LLkminus1 - LLk)
      k         <- length(centros)
      
      ## A     <- RAZAO * (1/(1-c)) * (N-k+1) * 1 * (phi/probMkminus1) ## ORIGINAL
      A     <- RAZAO * (1/(1-c)) * 1 * (phi/probMkminus1)
      
      alpha <- min(1, A)
      
      if(runif(1) < alpha){ 
        centros     <- new.centros
        particoes   <- new.particoes
        v_aceite[i] <- 1
      } 
    }
    
    ## Update step ------------------------------------------------------------
    if(passo == "Update"){
      ## Atualiza o vetor de médias
      medias    <- aggregate(y ~ particoes, FUN=mean)
      contagens <- as.data.frame(table(particoes))
      medias    <- merge(medias, contagens, by="particoes") 
      pesos     <- sigma2/( medias$Freq*sig2ma0 + sigma2 )
      
      ## Atualiza VETOR DE MEDIAS: "PRIORI" com "VEROSSIMILHANÇA LOCAL"
      vec.Medias[medias$particoes] <- rnorm(dim(medias)[1],
                                            mean = pesos*mu0 + (1-pesos)*medias$y,
                                            sd   = sqrt( 1/((1/sig2ma0) + (medias$Freq/sigma2)) )  )
      
      ## Atualiza o parâmetro de variância POR CONJUGAÇÃO...
      ## "PRIORI" com "VEROSSIMILHANÇA"
      ## Supondo médias conhecidas
      k <- length(centros)
      
      if(k < N){
        soma2  <- sum( (y - vec.Medias[particoes])^2 ) 
        ## Atualiza o parâmetro de variância POR MÁXIMA VEROSSIMILHANÇA...
        S2     <- soma2/(N-k)
        S2.new <- geoR::rinvchisq(n=1, df=N-k, scale=S2) ## PROBLEMA COM "df=N-k" quando "N-k -> 0"
        
        if(S2.new < 100){
          sigma2 <- S2.new
          v_aceite[i] <- 1
        }
      }
      v_aceite[i] <- 1
    }
    
    ## Shift step -------------------------------------------------------------
    if(passo == "Shift"){
      
      centros_shift_ok  <- centros
      new.centros       <- centros
      
      for(cntr in centros)    
      {
        ## Se todos os vizinhos de um centro também são centros: remove o centro do vetor
        if(all(viz[which(viz[,1]==cntr),2] %in% centros)) 
          centros_shift_ok <- centros_shift_ok[centros_shift_ok!=cntr]
      }
      
      ## Dentre os centros possíveis, escolhe um
      if(length(centros_shift_ok)==1){ 
        centro <- centros_shift_ok
      } else   centro <- sample(centros_shift_ok, size = 1)
      ## Vertor com os vizinhos do centro escolhido
      viz_do_centro <- viz[which(viz[,1]==centro),2]
      ## Vetor com os vizinhos livres do centro escolhido
      viz_livres     <- viz_do_centro[!(viz_do_centro %in% centros)]
      
      if(length(viz_livres)==1){ 
        shift <- viz_livres
      } else   shift <- sample(viz_livres, size = 1)
      
      new.centros[which(centros==centro)] <- shift
      
      n_Gk <- length(centros_shift_ok)
      m_Gj <- length(viz_livres)
      
      new.particoes  <- RcppPartition(viz, new.centros)
      
      vec.Medias[shift] <- vec.Medias[centro]
      
      new.centros_shift_ok <- new.centros
      ## Processo idêntico ao anterior para o novo vetor de centros
      for(nc in new.centros)   
      {
        if(all(viz[which(viz[,1]==nc),2] %in% new.centros)) 
          new.centros_shift_ok <- new.centros_shift_ok[new.centros_shift_ok!=nc]
      }
      
      ## Vizinhos do novo centro
      viz_new_centro  <- viz[which(viz[,1]==shift),2]
      ## Vizinhos livres do novo centro
      new.viz_livres  <- viz_new_centro[!(viz_new_centro %in% new.centros)]
      
      new.n_Gk <- length(new.centros_shift_ok)
      new.m_Gj <- length(new.viz_livres)
      
      ## Cálculo da razão de Verossimilhança
      LLk       <- sum( dnorm(y, mean = vec.Medias[particoes],     sd=sqrt(sigma2), log=TRUE ) )
      LLkshift  <- sum( dnorm(y, mean = vec.Medias[new.particoes], sd=sqrt(sigma2), log=TRUE ) )
      RAZAO     <- exp(LLkshift - LLk)
      
      A = RAZAO * (n_Gk/new.n_Gk) * (m_Gj/new.m_Gj)
      alpha <- min(1, A)
      
      if(runif(1) < alpha){
        centros     <- new.centros
        particoes   <- new.particoes
        v_aceite[i] <- 1
      }
    }
    
    ## Switch step ------------------------------------------------------------
    if(passo == "Switch"){
      
      ## Escolhe dois centros e troca
      switch      <- sample(1:length(centros), size=2)
      new.centros <- replace(centros, list = switch, centros[rev(switch)])
      
      new.particoes  <- RcppPartition(viz, new.centros)
      
      ## Cálculo da razão de Verossimilhança
      LLk       <- sum( dnorm(y, mean = vec.Medias[particoes],     sd=sqrt(sigma2), log=TRUE ) )
      LLkswitch <- sum( dnorm(y, mean = vec.Medias[new.particoes], sd=sqrt(sigma2), log=TRUE ) )
      RAZAO     <- exp(LLkswitch - LLk)
      
      A = RAZAO
      alpha <- min(1, A)
      
      if(runif(1) < alpha){
        centros     <- new.centros
        particoes   <- new.particoes
        v_aceite[i] <- 1
      }
      
    }
    
    ## - - - - - - - - - - - - - - - - - - - - -
    k_vector[i]  <- length(centros)
    v_sigma2[i]  <- sigma2
    v_passos[i]  <- passo
    mat.Yhat[,i] <- vec.Medias[particoes]
    v_centros[i] <- paste(centros, collapse=";")
    
    ## Set progress bar	
    setTxtProgressBar(pb, i, "running Reversible Jump...")
    
    ## Preencher matriz de frequência
    if(i > burn_in) freq_matrix <- freq_matrix + RcppFreqMatrix(particoes)
    
  }
  
  Sys.sleep(1)
  close(pb)  # Ends progress bar
  
  ## Results processing -------------------------------------------------------
  
  # Remove burn-in
  seq.burn <- -(1:burn_in)
  k_vector <- k_vector[seq.burn]
  v_sigma2 <- v_sigma2[seq.burn]
  v_passos <- v_passos[seq.burn]
  mat.Yhat <- mat.Yhat[,seq.burn]
  v_aceite <- v_aceite[seq.burn]
  v_centros <- v_centros[seq.burn]
  
  ## Médias "a posteriori" e "intervalos de credibilidade"
  Yhat <- apply(mat.Yhat, MARGIN=1, FUN=median)
  lwr  <- apply(mat.Yhat, MARGIN=1, FUN=function(x) quantile(x, probs=0.05))
  upr  <- apply(mat.Yhat, MARGIN=1, FUN=function(x) quantile(x, probs=0.95))
  
  ## Processamento da Matrix de Frequências de "conexões"
  ## Pegando o resultado do Leandro
  matConexoes <- freq_matrix
  
  ## Preenche TODA a matriz
  matConexoes[lower.tri(matConexoes)] <- t(matConexoes)[lower.tri(matConexoes)]
  maximo <- max(as.vector(matConexoes), na.rm=TRUE)
  matConexoes       <- maximo - matConexoes
  diag(matConexoes) <- maximo + 1
  
  ## Forma as partições utilizando análise Hierárquica de clusters
  clusters <- hclust(as.dist(matConexoes), c("single","complete")[1])
  
  output <- list(mean.info    = cbind(lwr, Yhat, upr),
                cluster.info = clusters,
                matConexoes  = matConexoes,
                k.MCMC       = k_vector,
                mean.y       = mean.y,
                sd.y         = sd.y,
                vec.centros  = v_centros)
  
  return(output)
}
