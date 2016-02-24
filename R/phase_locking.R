#' pahseLocking calculation for two vectors
#'
#' pahseLocking is written for testing purposes.
#' if we are calculating phase locking for a large matrix,
#' running it per pair is slow
#' @param x - a vector, should be a time series vector
#' @param y - a vector, should be a time series vector
#' @param n - time period tatio
#' @param m - time period tatio
#' @return  list with 1. entropy rho 2. gamma  3. Strobo index lambda 4. Strobo index angle
#' @examples
#' x <- sin(1:200)
#' y <- sin(2*1:200)
#' n=1
#' m=2
#' phaseLocking(x, y, n=1, m=2)
#' @references
#' [1]
#' @export
phaseLocking <- function(x, y, n = 1, m = 1) {
  if (length(x) != length(y))
    stop(paste("Expression length is different",
               "please check the data"))
  
  # Hilbert transformation from RSEIS
  x_hilb = RSEIS::hilbert(RSEIS::detrend(x))
  y_hilb = RSEIS::hilbert(RSEIS::detrend(y))
  
  # calculate instantaneous phase
  x_phase = RSEIS::LocalUnwrap(Arg(x_hilb), cutoff = pi)
  y_phase = RSEIS::LocalUnwrap(Arg(y_hilb), cutoff = pi)
  
  p2p = 2 * pi
  x_mod = x_phase %% p2p
  y_mod = y_phase %% p2p
  
  ##n:m locking
  lock = n * x_phase - m * y_phase
  
  # We test three measures of syncrhonization
  #  1. entropy rho
  #  2. gamma
  #  3. Strobo index lambda
  
  #1. entropy rho
  #   index according to Shannon's entropy
  
  #transfer to [0 2*pi)
  lock_mod = lock %% p2p
  
  #transfer to [-pi pi)'
  lock_mod [lock_mod >= pi] = lock_mod [lock_mod >= pi]  - p2p
  
  #how is the nbins (number of bins) calculated?
  nbins = round(exp(0.626 + 0.4 * log(length(x) - 1)))
  hist_bins = (0:nbins) * p2p / nbins - p2p / 2
  
  #bin the phase
  bin_result = .bincode(lock_mod,hist_bins)
  bin_result = bin_result[bin_result != 0]
  bin_fraction = prop.table(bin_result)
  
  SS = -sum(bin_fraction * log(bin_fraction))
  Smax = log(nbins)
  rho = (Smax - SS) / Smax
  
  # 2. gamma
  #    intensity of the first Fourier mode of the distrbution
  gamma = (mean(cos(lock))) ^ 2 + (mean(sin(lock))) ^ 2  # which is lmb^2;
  
  # 3. strobo index lambda
  #    not appropriate for small # of datea points. %%
  
  nn = hist(x_mod, breaks = (0:nbins) * p2p / nbins, plot = FALSE)$counts  ###Need to be start from 0
  nb = length(nn[nn > 0])
  bin_index = nn > 0
  
  lambda1 = rep(0,nbins)
  lambda2 = rep(0,nbins)
  
  sortObj = sort(x_mod, index.return = TRUE)
  pd1s = sortObj$x
  idx = sortObj$ix
  
  counter1 = 0;
  for (ii in 1:nbins) {
    theta1 = (ii - 1) * p2p / nbins + p2p / 2 / nbins
    for (jj in 1:nn[ii]) {
      counter1 = counter1 + 1;
      lambda1[ii] = lambda1[ii] + exp(complex(imaginary = y_mod[idx[counter1]]))
      lambda2[ii] = lambda2[ii] + exp(complex(imaginary = y_mod[idx[counter1]] -
                                                x_mod[idx[counter1]]))
    }
    lambda1[ii] = lambda1[ii] / nn[ii]
    lambda2[ii] = lambda2[ii] / nn[ii]
  }
  
  lmb1 = sum(abs(lambda1[bin_index])) / nb
  lmb2 = sum(abs(lambda2[bin_index])) / nb
  ang2 = mean(Arg(lambda2[bin_index]));
  #  #4, a simple index by the mean phase difference.
  # there is a problem, the phase difference between two random time series
  # follows roughly a normal distribution, the circular mean of the phase
  # difference is then biased at 0, with significant amplitude (~0.67)
  lam = mean(exp(complex(
    real = rep(0, length(y)), imaginary = lock
  )));
  lmb = abs(lam);
  ang = Arg(lam);
  
  return(list(
    "entropy rho" = rho,
    "gamma" = gamma,
    "strobo lmb" = (lmb1 + lmb2) / 2,
    "strobo ang" = ang
  ))
}

#' pahseLocking calculation for matrix
#'
#' pahseLocking is written for testing purposes.
#' if we are calculating phase locking for a large matrix,
#' running it per pair is slow
#' @param x - a matrix, rows are genes, columns are sorted by time
#' @param n - time period tatio
#' @param m - time period tatio
#' @param method - one in c(entrophy, gamma, strobo)
#' @return   list with 4 matrices
#' @return 1. entropy rho 2. gamma  3. Strobo index lambda 4. Strobo index angle
#' @examples 
#' x <- replicate(100, rnorm(20))
#' n=1
#' m=1
#' phaseLockingMatrix(x, n, m)
#' @references
#' [1]
#' @export

phaseLockingMatrix <- function(x, n = 1, m = 1,method = "strobo") {
  #remove rows with NA
  x = x[complete.cases(x),]
  
  rhoM = matrix(0,nrow = nrow(x),ncol = nrow(x))
  gammaM = matrix(0,nrow = nrow(x),ncol = nrow(x))
  strobolmb = matrix(0,nrow = nrow(x),ncol = nrow(x))
  stroboang = matrix(0,nrow = nrow(x),ncol = nrow(x))
  

  
  pb <- txtProgressBar(min = 0, max = nrow(x) - 2, style = 3)
  for (i in 1:(nrow(x) - 2)) {
    setTxtProgressBar(pb, i)
    this_series = x[i,]
    
    phase = apply(x[1:nrow(x),], 1, function(y)
      phaseLocking(this_series,y,n,m))
    
    rhoM[i,] = unlist(lapply(phase, `[[`, 1))
    gammaM[i,] = unlist(lapply(phase, `[[`, 2))
    strobolmb[i,] = unlist(lapply(phase, `[[`, 3))
    stroboang[i,] = unlist(lapply(phase, `[[`, 4))
  }
  close(pb)
  
  #last pair
  i = i + 1
  phase = phaseLocking(x[i,],x[i + 1,],n,m)
  rhoM[i,i + 1] = phase[[1]]
  gammaM[i,i + 1] = phase[[2]]
  strobolmb[i,i + 1] = phase[[3]]
  stroboang[i,i + 1] = phase[[4]]
  
  diag(rhoM) <- 1
  diag(gammaM) <- 1
  diag(strobolmb) <- 1
  diag(stroboang) <- 1
  
  #make the NA's to be 0
  rhoM[is.na(rhoM)] <- 0
  gammaM[is.na(gammaM)] <- 0
  strobolmb[is.na(strobolmb)] <- 0
  stroboang[is.na(stroboang)] <- 0
  
  #Take the maximum value, make the matrix symetric
  temp1 = t(rhoM)
  ind = which(rhoM < temp1)
  rhoM[ind] = temp1[ind]
  
  temp2 = t(gammaM)
  ind = which(gammaM < temp2)
  gammaM[ind] = temp2[ind]
  
  temp3 = t(strobolmb)
  ind = which(strobolmb < temp3)
  strobolmb[ind] = temp3[ind]
  
  #
  rownames(rhoM) <- colnames(rhoM) <- rownames(gammaM) <-
    colnames(gammaM) <- rownames(strobolmb) <- colnames(strobolmb) <-
    rownames(stroboang) <- colnames(stroboang) <- rownames(x)
  
  #Add a small number to avoid 0 for
#   return(
#     list(
#       "entropy_rho" = rhoM + 10 ^ (-6),
#       "gamma" = gammaM + 10 ^ (-6),
#       "strobo_lmbda" = strobolmb + 10 ^ (-6),
#       "strobo_angle" = stroboang
#     )
#   )
  
  return(
    list(
      "entropy_rho" = rhoM,
      "gamma" = gammaM ,
      "strobo_lmbda" = strobolmb,
      "strobo_angle" = stroboang
    )
  )
}

#' Granger cause calculation for two vectors
#' @param x - a matrix, rows are genes, columns are times
#' @param p - lag
#' @return  significant p-value with hypothesis that x granger cause y with lag=order
#' @examples
#' #example 1
#' x = replicate(100, rnorm(5))
#' rownames(x)=paste("gene",1:nrow(x))
#' grangerTest(x, p=10)
#'
#' #example 2
#' x = rbind(sin(1:100),cos(1:100))
#' rownames(x)=c("sin","cos")
#' grangerTest(x, p=4)
#'
#' @references
#' [1] http://www.r-bloggers.com/chicken-or-the-egg-granger-causality-for-the-masses/
#' [2] Thurman W.N. & Fisher M.E. (1988), Chickens, Eggs, and Causality, or Which Came First?, American Journal of Agricultural Economics, 237-238.
#' @export

grangerTest = function(x,p = 1) {
  #remove rows with NA
  x = x[complete.cases(x),]
  
  #detrend for the signal
  x = t(apply(x,1,function(x)
    detrend(x)))
  
  a = MSBVAR::granger.test(t(x),p = p)
  F_statistics = a[,1]
  p_value = a[,2]
  n = nrow(x)
  
  #assign NAs to be 0 so to avoid error in next matrix forming step
  F_statistics[is.na(F_statistics)]<-0
  p_value[is.na(p_value)]<-1
  
  #F_statistics and p_value are vectors, need to change to matrix
  #add some NA at the diagonal places
  #store the index you want to insect 0 or 1
  na_ind=seq(1, n^2, by = n+1) 
  temp_all <- numeric(length(F_statistics)+length(na_ind)) 
  temp_all[na_ind] <- NA 
  
  temp_all[!is.na(temp_all)] <- F_statistics 
  F_statistics_matrix=matrix(temp_all,nrow=n,ncol=n)
  #Since we are not considering self loop, any value works
  F_statistics_matrix[is.na(F_statistics_matrix)]<-0
  rownames(F_statistics_matrix) = rownames(x)
  colnames(F_statistics_matrix) = rownames(x)
  
  temp_all[!is.na(temp_all)] <- p_value
  p_value_matrix=matrix(temp_all,nrow=n,ncol=n)
  #Since we are not considering self loop, any value works
  p_value_matrix[is.na(p_value_matrix)]<-0
  rownames(p_value_matrix) = rownames(x)
  colnames(p_value_matrix) = rownames(x)
  
  
  #Take the maximum value, make the matrix symetric
  temp1 = t(F_statistics_matrix)
  ind = which(F_statistics_matrix < temp1)
  F_statistics_matrix[ind] = temp1[ind]
  
  #Take the minimum value, make the matrix symetric
  temp2 = t(p_value_matrix)
  ind = which(p_value_matrix > temp2)
  p_value_matrix[ind] = temp2[ind]
  
  return(list(
    "F_statistics_matrix" = F_statistics_matrix,"p_value_matrix" = p_value_matrix
  ))
}


#' Coherence for time series matrix
#' @param x - a matrix, rows are genes, columns are times
#' @return  maximun coherence matrix between each pair of time series
#' @examples
#' #example 1
#' x = replicate(100, rnorm(5))
#' coherenceTest(x)
#'
#' #example 2
#' x = rbind(sin(1:100),cos(1:100))
#' rownames(x)=c("sin","cos")
#' coherenceTest(x)
#'
#' @references
#'
#'
#' @export
coherenceTest <- function(x)
{
  #remove rows with NA
  x = x[complete.cases(x),]
  
  #detrend for the signal
  x = t(apply(x,1,function(x)
    detrend(x)))
  
  result = matrix(0,nrow = nrow(x),ncol = nrow(x))
  
  if (nrow(x) < 2)
    
    return("You need to have at least 2 time series to run coherence test")
  
  if (nrow(x) == 2)
    return(max(
      spec.pgram(
        cbind(x[1,], x[2,]), fast = FALSE, taper = FALSE,
        spans = c(3,3),plot = FALSE
      )$coh
    ))
  pb <- txtProgressBar(min = 0, max = nrow(x) - 2, style = 3)
  for (i in 1:(nrow(x) - 2))
  {
    # update progress bar
    setTxtProgressBar(pb, i)
    this_series = x[i,]
    
    #return of coh(): 1st column - frequency, 2nd column - coherence
    #Since cohrence is symmeteric, we only calculate the upper triangle
    coherence = apply(x[(i + 1):nrow(x),], 1, function(y)
      spec.pgram(
        cbind(this_series, y), fast = FALSE, taper = FALSE,
        spans = c(3,3),plot = FALSE
      )$coh)
    result[i,(i + 1):nrow(x)] = apply(coherence, 2, max)
  }
  
  
  close(pb)
  
  #the last pair
  result[i + 1,i + 2] = max(spec.pgram(
    cbind(x[i + 1,], x[i + 2,]), fast = FALSE, taper = FALSE,
    spans = c(3,3),plot = FALSE
  )$coh)
  
  #assign NaN to be 0
  result[is.nan(result)] <- 0
  
  #some elements may return a number bigger than 1 due to machine error
  #make them to be 1
  result[result > 1] <- 1
  
  #make the result matrix symmetric
  result[lower.tri(result)] = t(result)[lower.tri(result)]
  diag(result) <- 1
  
  rownames(result) = rownames(x)
  colnames(result) = rownames(x)
  return("coh" = result)
  
}
