#' GWAS with Andersen-Gill score test
#'
#' @param x a bedmatrix
#' @param null.model an object of class nullModel
#' @param beg Index of the first SNP tested for association
#' @param end Index of the last SNP tested for association
#'
#' @details 
#' 
#' @return A data frame giving for each SNP the association statistics.
#' 
#' @examples # TBC
#' 
#' @export
andersen.gill.score_TEST <- function(x, null.model, beg = 1, end = ncol(x)) {

  if(beg < 1 | end > ncol(x)) stop("range too wide")
  if(is.null(x@mu) | is.null(x@p)) stop("Need mu and p to be set in x (use set.stats)")

  # check ids ...
  id <- as.factor(null.model@ids)
  times <- rle(null.model@ids)$length  # number of lines for each id
  if(nlevels(id) != length(times)) {
    stop("The null model should be constructed with all events related to an individual in contiguous lines")
  }
  if(any(unique(null.model@ids) != x@ped$id)) {
    stop("The null model should be constructed with the indivuals in the same order than in the bed matrix")
  }
    
  if(null.model@float) 
    t <- GWAS_score_test_AG_bed_float_TEST(x@bed, x@p, times, beg - 1L, end - 1L, null.model@ptr)
  else
    t <- GWAS_score_test_AG_bed_double_TEST(x@bed, x@p, times, beg - 1L, end - 1L, null.model@ptr)

  t$p <- pchisq(t$score, df = 1, lower.tail = FALSE)

  # mise en forme
  L <- data.frame(chr = x@snps$chr, pos = x@snps$pos, id = x@snps$id, A1 = x@snps$A1, A2 = x@snps$A2, freqA2 = x@p, stringsAsFactors = FALSE)
  if(beg > 1 | end < ncol(x))  # avoid copy
    L <- L[beg:end,] 

  data.frame( c( L, t) )
}

