
power_z <- function(r, n, c = 0, 
                    compare = FALSE, 
                    alpha = 0.05){
  
  # fisher z
  z <- atanh(abs(r))
  
  # variance of z
  var_z <- 1/ (n - c - 3)
  
  # differnece ?
  if(compare == TRUE){
    
    var_z <- var_z * 2
    
  }
  
  # z score
  z_score <- z/ sqrt(var_z)
  
  # quantile
  q <- stats::qnorm(1 - alpha / 2)
  
  # power
  1 - stats::pnorm(q - z_score)
  
}


power_eq <- function(r_lb = 0.1, r_ub = 0.1, n, k,
                     compare = TRUE, alpha = 0.05){
  
  
  if(compare == TRUE){
    
    z_var <- (1 / (n - k - 3)) * 2
    
    } else {
    
    z_var <- (1 / (n - k - 3)) 
    
    }
  
  
  z_lb <- abs(atanh(r_lb)) / sqrt(z_var)
  z_ub <- abs(atanh(r_ub)) / sqrt(z_var)
  
  pwr_lb <- 2 * (pnorm(z_lb - qnorm(1-alpha)) + pnorm(-z_lb - qnorm(1-alpha))) - 1
  pwr_ub <- 2 * (pnorm(z_ub - qnorm(1-alpha)) + pnorm(-z_ub - qnorm(1-alpha))) - 1
  
  min(pwr_lb, pwr_ub)
  
}


ordinal_power <- function(n, pcor_mat, 
                          levels = 5, 
                          probs = NULL, 
                          emp = TRUE, B = 500){
  
  
  adj <- ifelse(pcor_mat == 0, 0, 1)
  
  which_nonzero <- which(adj[upper.tri(adj)] == 1)
  
  temp <- pcor_mat * -1
  
  diag(temp) <- 1
  
  cor_mat <- cov2cor( solve(temp)) 
  
  Y <- replicate(50, gen_ordinal_new(n = n, 
                                     cor_mat = cor_mat, 
                                     levels = levels, 
                                     probs = probs, emp = TRUE))
  
  # cluster
  cl <- parallel::makeCluster(4) 
  
  # register cluster
  doParallel::registerDoParallel(cl)
  
  # cluster export
  # parallel::clusterExport(cl, varlist  = c("n", "Y", "cor_mat"))
  
  mse <- parallel::parLapply(cl = cl, 1:25, function(x) 
  {Y_temp <-  as.data.frame( Y[,,x] )
  colnames(Y_temp) <- 1:ncol(Y_temp)
  cors <-  lavaan::lavCor(Y_temp);
  class(cors) <- "matrix"
  mean((cors[upper.tri(cors)] - cor_mat[upper.tri(cors)])^2)
  }
  )
  # data closest to the true matrix
  dat <- as.data.frame( Y[,,which.min(mse)] )
  
  colnames(dat) <- 1:ncol(pcor_mat)
  # cluster
  cl <- parallel::makeCluster(4) 
  
  # register cluster
  doParallel::registerDoParallel(cl)
  
  
  
  test <- parallel::parLapply(cl, 1:B, function(x)  
  {cors <- lavaan::lavCor(dat[sample(1:n, size = n, replace = T), ], 
                          ordered = names(dat));
  class(cors) <- "matrix"
  pcors <- cov2cor( solve(cors)) * -1;
  atanh(pcors[upper.tri(pcors)] );
  
  } )
  
  
  z_se <- apply(na.omit(do.call(rbind, test)), 2, sd)[which_nonzero]
  
  pcors <- pcor_mat[upper.tri(pcor_mat)][which_nonzero]
  
  pwr_s <- 1 - stats::pnorm(1.96 - (abs(atanh(pcors)) / z_se))
  
  pwr <- 1 - stats::pnorm(1.96 - (abs(atanh(pcors)) / mean(z_se)))
  
  return_object<- list(pwr = pwr, 
                       pwr_s = pwr_s, z_se = z_se, 
                       non_zero = which_nonzero, 
                       dat = dat, mse = mse)
  
}
















gen_ordinal_new <- function(n,  cor_mat, 
                            levels, 
                            probs = NULL,  
                            emp = FALSE){
  
  p <- ncol(cor_mat)
  
  ls <- list()
  
  if(is.null(probs)){
    
    probs <- rep(1, levels)
    
  } else{
    
    if(length(probs) != levels){
      stop("length of probs and levels must be the same")
    }
    
    
  }
  
  for(i in 1:p){
    
    temp <- table(sample(c(1:levels), 
                         size = n, 
                         replace = T, prob = probs)) 
    ls[[i]] <- as.numeric(temp / sum(temp))
  }
  
  junk <- capture.output( data <- rmvord_naiv(n = n,
                                              probs =  ls,
                                              Cor = cor_mat, 
                                              emp = emp))
  
  
  data
  
}


rmvord_naiv <- function(n, probs, Cors, emp) {
  q = length(probs)
  categ_probs = 0
  cumul_probs = list(0)
  quant_probs = list(0)
  means = 0
  vars = 0
  var.wt = function(x, w) {
    m = weighted.mean(x = x, w = w)
    sum((x[1:length(x)] - m)^2 * w[1:length(x)])
  }
  for (i in 1:q) {
    categ_probs[i] = length(probs[[i]])
    cumul_probs[[i]] = cumsum(1:categ_probs[i]/10^12 + probs[[i]])
    cumul_probs[[i]][categ_probs[i]] = 1
    quant_probs[[i]] = qnorm(p = cumul_probs[[i]], mean = 0, 
                             sd = 1)
  }
  retval = MASS::mvrnorm(n = n, 
                         mu  = rep(0,q), 
                         Sigma =  Cors, 
                         empirical = emp)
  for (i in 1:q) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quant_probs[[i]]), 
                      right = FALSE)
  }
  retval
}

