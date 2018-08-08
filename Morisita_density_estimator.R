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
  
  if (n < 3)
    stop("n must be equal to or larger than 3.")
  N <- nrow(r) # number of the sampling points
  k <- ncol(r) # number of the directions
  
  K <- sum(!is.na(r))   # total number of counts excluding NA
  k_prime <- apply(!is.na(r), 1, sum)
  m_hat1 <- k * (n - 1) / K * sum(1 / r^2, na.rm = TRUE)
  m_hat2 <- sapply(1:N,
                   function(i)
                     k * (n * k_prime[i] - 1) / sum(r[i, ]^2, na.rm = TRUE))
  m_hat2 <- sum(m_hat2) / N
  m_hat0 <- (m_hat1 + m_hat2) / 2
  if (m_hat1 > m_hat2)
    return(m_hat1 / pi)
  else
    return(m_hat0 / pi)
}
