# Morisita's Angle-order density estimator (Morisita 1957)
# Morisita M. (1957) A new method for the estimation of density by spacing method
# applicable to nonrandomly distributed populations. Physiology and Ecology 7:134-144
# (in Japanese).
# http://reference.morisita.or.jp/paper_pdf/29.pdf
# English translation: Forest Service translation number 11116,
# USDA Forest Service, Washington, D.C., USA.

estimate_density <- function(r, n = 3L) {
  # r: a matrix containing measurement data (distance to n-th neibouring point)
  #    The rows correspond to the sampling points, and the columns correspond to
  #    the directions.
  # n: an integer indicating that n-th neighbour is measured (n >= 3)
  # Returns the estimated density of the population
  
  if (!is.matrix(r))
    stop("r must be a matrix.")
  if (n < 3)
    stop("n must be equal to or larger than 3.")
  N <- nrow(r) # number of the sampling points
  if (N < 2)
    stop("The number of the sampling points must be equal to or larger than 2.")
  k <- ncol(r) # number of the directions
  if (k < 2)
    stop("The number of the directions must be equal to or larger than 2.")
  
  K <- sum(!is.na(r))   # total number of measurements excluding NA
  k_prime <- rowSums(!is.na(r)) # number of measurements for each sampling point
  m_hat1 <- k * (n - 1) / K * sum(1 / r^2, na.rm = TRUE)
  m_hat2i <- k * (k_prime * n - 1) / rowSums(r^2, na.rm = TRUE)
  m_hat2 <- sum(k_prime * m_hat2i) / K
  m_hat0 <- (m_hat1 + m_hat2) / 2
  density <- ifelse(m_hat1 > m_hat2, m_hat1 / pi, m_hat0 / pi)

  return(density)
}
