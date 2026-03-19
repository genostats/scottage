require(scottage)
source("my.null.model.r")
source("my.score.test.r")

data <- readRDS("mydata.rds")
formula <- Surv(start, stop, etat) ~ cluster(id)

# -------------------------------------------

nm0 <- my.null.model(formula, data)
s0 <- my.score.test(nm0, data$tt)

scottage(use.float = FALSE)
nm1 <- null.model(formula, data)
s1 <- score.test(nm1, data$tt)

s2 <- coxph(Surv(start, stop, etat) ~ tt + cluster(id), data, ties = "breslow")$rscore

stopifnot( abs(s0 - s1) < 1e-6 )
stopifnot( abs(s0 - s2) < 1e-6 )
stopifnot( abs(s1 - s2) < 1e-6 )

scottage(use.float = TRUE)
nm3 <- null.model(formula, data)
s3 <- score.test(nm3, data$tt)
stopifnot( abs(s0 - s1) < 1e-5 )

# ---- un test avec des valeurs manquantes dans tt ----

X <- data$tt
X[1:5] <- NA
s0.na <- my.score.test(nm0, X)
s1.na <- score.test(nm1, X)
stopifnot( abs(s0.na - s1.na) < 1e-6 )

# -------------------------------------------

formula2 <- Surv(start, stop, etat) ~ covar + covar2 + cluster(id) 

nm4 <- my.null.model(formula2, data)
s4 <- my.score.test(nm4, data$tt)

scottage(use.float = FALSE)
nm5 <- null.model(formula2, data)
s5 <- score.test(nm5, data$tt)

stopifnot( abs(s5 - s4) < 1e-6 )

