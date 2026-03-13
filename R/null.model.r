#' Creates null model object
#'
#' @param formula a formula, whose response is a survival object (a call to the 'Surv' function)
#' @param data a data frame
#'
#' @export
null.model <- function(formula, data) {
  # fit du modele
  fit <- survival::coxph(formula, ties = "breslow", data = data)
  
  # récuperation de start, stop, status
  start  <- as.numeric(fit$y[, "start"])
  stop   <- as.numeric(fit$y[, "stop"])
  status <- as.numeric(fit$y[, "status"])
  
  # récupération de la variable cluster
  cl <- as.character(fit$call$cluster)
  if(!is.null(cl) && length(cl) > 0) {
    id <- data[, cl]
  } else { # pas de cluster : un id unique par ligne
    id <- 1:nrow(data)   
  }
  id <- as.factor(id)
  
  
  Terms <- terms(formula, specials = "cluster", data = data)
  tcl <- attr(Terms, "specials")$cluster
  
  # récupérer les covariables normales (excluant response et cluster )
  newTerms <- delete.response(Terms)
  term.labels <- attr(newTerms, "term.labels")
  normal.terms <- term.labels[!grepl("^cluster", term.labels)]
  
  if(length(normal.terms) > 0) {
    new.formula <- as.formula(paste("~", paste(normal.terms, collapse = "+")))
    # préparation de la base pour le modèle prenant en compte les calculs comme log, sqrt ..
    mf <- model.frame(new.formula, data)
    Z <- model.matrix(new.formula, mf)
    ## A quel terme de la formule appartient chaque colonne de la matrice de design
    attr(Z, "assign") 
    ## récupérer tous les attributs de Z
    Xatt <- attributes(Z)
    ## regrouper les colonnes de la matrice de design Z par terme de la formule
    assign <- attrassign(Z, Terms) 
    ## supression de l'intercept
    adrop <- 0
    xdrop <- Xatt$assign %in% adrop
    Z <- Z[, !xdrop, drop = FALSE]
    no.covars <- FALSE
  } else {
    Z <- matrix(nrow = nrow(data), ncol = 0)
    no.covars <- TRUE
  }
  
  # récupération de beta
  beta <- fit$coefficients
  
  # temps d'événement 
  T <- sort( stop[ status == 1] )
  T <- T[ !duplicated(T) ]
  
  # pour chaque ligne du tableau, les indices des temps d'événement
  S <- match(stop, T)
  
  m <- length(T)
  n <- length(start)
  p <- ncol(Z)
 
  # Pré calcul des liste d'indices
  # indices des événements
  events <- which(status == 1)
  # indices des temps qui sont dans l'intervalle (start, stop], pour chaque ligne
  L.times <- vector("list", n)
  for(jk in 1:n) {
    L.times[[jk]] <- which(start[jk] < T & T <= stop[jk])
  }
  L.events <- vector("list", n)
  # indices des lignes avec un événement dans l'intervalle (start, stop]
   for(jk in 1:n) {
    L.events[[jk]] <- which(status == 1 & start[jk] < stop & stop <= stop[jk])
  }

  # On calcule en décomposant chaque Yj (Zj) en une somme de Y_jk (Z_jk) qui correspondent
  # chacun à une ligne du tableau
  
  if(no.covars) { # on calcule juste S0
    S0 <- numeric(m)
    for(jk in 1:n) {
      S0[ L.times[[jk]] ] <- S0[ L.times[[jk]] ] + 1
    }
    # et c'est fini pour ce cas
    return(nullModel(formula = formula, data = data, start = start, stop = stop, 
                     status = status, id = id, ev.times = T, ev.index = S, 
                     S0 = S0, exp.beta.Z = matrix(1, ncol = 1, nrow = n), 
                     W = matrix(0, nrow = 0, ncol = length(unique(id))), 
                     inverse.variance = matrix(0, nrow = 0, ncol = 0), 
                     events = events, L.times = L.times, L.events = L.events))
  }
  
  exp.beta.Z <- exp(as.numeric(Z %*% beta))
  Z.exp.beta.Z <- Z * exp.beta.Z
  
  S0 <- numeric(m)
  S1 <- matrix(0, nrow = p, ncol = m) # on met les dérivées sucessives en lignes
  for(jk in 1:n) {
    S0[ L.times[[jk]] ] <- S0[ L.times[[jk]] ] + exp.beta.Z[jk]
    S1[, L.times[[jk]] ] <- S1[, L.times[[jk]] ] + Z.exp.beta.Z[jk,]
  }
  
  # -------------------------------------------------
  # score par rapport aux covariables
  # on calcule la contribution au score de la ligne jk ... 
  U <- matrix(0, nrow = p, ncol = n)
  for(jk in events) {
    s <- S[jk]
    U[,jk] <- Z[jk,] - S1[,s]/S0[s]
  }
  
  # et EU, de même, ligne à ligne
  EU <- matrix(0, nrow = p, ncol = n)
  for(jk in 1:n) {
    for(l in L.events[[jk]]) {
      s <- S[l]
      EU[,jk] <- EU[,jk] + exp.beta.Z[jk]/S0[s] * (Z[jk,] - S1[,s]/S0[s])
    }
  }
  
  # Attention c'est toujours ligne à ligne
  W <- U - EU
  
  # On regroupe pour avoir les valeurs par individus
  tmp <- tapply(1:n, id, \(I) rowSums(W[, I, drop = FALSE]) )
  W <- matrix(unlist(tmp), nrow = p)
  
  # et voila la matrice de variance robuste
  Vkk <- tcrossprod(W)
  
  nullModel(formula = formula, data = data, start = start, stop = stop, 
            status = status, id = id, ev.times = T, ev.index = S, 
            S0 = S0, exp.beta.Z = exp.beta.Z, W = W, inverse.variance = solve(Vkk),
            events = events, L.times = L.times, L.events = L.events)
}

