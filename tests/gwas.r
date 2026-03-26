require(scottage)
data <- readRDS("mydata.rds")
formula <- Surv(start, stop, etat) ~ cluster(id)

# on fabrique une bed matrix de 1 colonne...
G <- data$tt[ !duplicated(data$id) ]
u.id <- unique(data$id)
x <- as.bed.matrix( matrix(G, ncol = 1), 
                    fam = data.frame(famid = u.id, id = u.id, father = 0, mother = 0, sex = 0, pheno = 0),
                    bim = data.frame(chr = 1, id = "rs01", dist = NA, pos = 10, A1 = "A", A2 = "C"))


# --------- test en double -------------
scottage(use.float = FALSE)
# avec scottage::null.model
nm1 <- null.model(formula, data)
s1 <- score.test(nm1, data$tt)

# et avec la bed matrix
s2 <- andersen.gill.score(x, nm1)

stopifnot( abs(s1 - s2$score) < 1e-15 )


# --------- test en float -------------
scottage(use.float = TRUE)
nm3 <- null.model(formula, data)
s3 <- score.test(nm3, data$tt)

# et avec la bed matrix
s4 <- andersen.gill.score(x, nm3)

stopifnot( abs(s3 - s4$score) < 1e-15 )



