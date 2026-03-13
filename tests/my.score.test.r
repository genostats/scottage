my.score.test <- function(null.model, X) {
  start  <- null.model$start
  stop   <- null.model$stop
  status <- null.model$status
  id <- null.model$id

  T <- null.model$ev.times
  S <- null.model$ev.index

  S0 <- null.model$S0

  exp.beta.Z <- null.model$exp.beta.Z
  W <- null.model$W
  Vkk_i <- null.model$inverse.variance

  m <- length(T)
  n <- length(start)

  if(length(X) != n) stop("length(X) is not ok")

  # -------------------------------------------------------------------------------------
  # Calcul de la composante du score p/r à X
  # On a besoin de la dérivée p/r à gamma [le coefficient de X, qui est pris à 0
  # dans tous les calculs]

  S1.gamma <- numeric(m)
  for(jk in 1:n) {
    interval <- (start[jk] < T & T <= stop[jk])
    S1.gamma[interval] <- S1.gamma[interval] + X[jk] * exp.beta.Z[jk] ## diff avec s1
  }

  # on calcule la contribution au score de la ligne jk ... 
  # vectorisable mais comme tout ça va passer en C++
  Ug <- numeric(n)
  for(jk in 1:n) {
    if(status[jk] == 1) { # ou encore !is.na(S[jk]) ... 
      s <- S[jk]
      Ug[jk] <- X[jk] - S1.gamma[s]/S0[s]
    }
  }

  # et EU, de même, ligne à ligne
  EUg <- numeric(n)
  for(jk in 1:n) {
    sta <- start[jk]
    sto <- stop[jk]
    for(l in 1:n) {
      if(status[l] == 0 || stop[l] <= sta || stop[l] > sto) next
      s <- S[l]
      EUg[jk] <- EUg[jk] + exp.beta.Z[jk]/S0[s] * (X[jk] - S1.gamma[s]/S0[s])
    }
  }

  # le vecteur W, toujours ligne à ligne
  Wg <- Ug - EUg

  # on passe à la valeur par individus
  Wg <- as.vector(tapply(Wg, id, sum))

  # Et enfin voilà la composante de Vgg de la variance
  Vgg <- sum(Wg**2)

  # Et Vkg [la covariance entre les deux morceaux]
  Vkg <- W %*% Wg

  # le score
  score <- sum(Wg)**2 / (Vgg - t(Vkg) %*% Vkk_i %*% Vkg)

  c("chi2" = score)
}

