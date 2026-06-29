#' Circular k-means clustering
#'
#' k-means on the torus \eqn{T^d}: a Lloyd iteration whose dissimilarity is the
#' summed cosine distance \eqn{\sum_j \{1 - \cos(\theta_j - \mu_j)\}} and whose
#' centres are the per-coordinate circular means
#' \eqn{\mu_j = \mathrm{atan2}(\overline{\sin\theta_j}, \overline{\cos\theta_j})}.
#' On angular data this is the right analogue of \code{\link[stats]{kmeans}}: it
#' respects wrap-around (\eqn{\theta} and \eqn{\theta + 2\pi} are the same point),
#' which Euclidean k-means on the raw radians does not, and its centres stay
#' \emph{on} the circle (unit resultant) rather than at the radially shrunk
#' arithmetic mean of the \eqn{(\cos\theta, \sin\theta)} embedding, so points are
#' assigned by true angular distance. For \eqn{d > 1} columns the distance sums
#' over coordinates, clustering on the product of circles -- the same torus
#' factorisation \code{\link{circ_mix}} uses for a joint angular response, which is
#' why \code{circ_mix} seeds and splits its components with this routine.
#'
#' Each Lloyd update sets a centre to the circular mean of its members, which is
#' exactly the minimiser of that cluster's summed cosine distance -- the circular
#' mean direction \eqn{\phi = \mathrm{atan2}(\sum\sin\theta_i, \sum\cos\theta_i)}
#' maximises \eqn{\sum_i \cos(\theta_i - \mu)}. So the alternation is coordinate
#' descent on one objective (monotone, convergent), the circular counterpart of
#' "the arithmetic mean minimises squared-Euclidean distance" behind ordinary
#' k-means. Equivalently, since \eqn{1 - \cos(\theta - \mu) = \tfrac12 \|e^{i\theta}
#' - e^{i\mu}\|^2}, this is spherical k-means on the unit circle with centres
#' projected back onto it -- the hard-assignment limit of a von Mises mixture with
#' common concentration (Banerjee et al., 2005), which is why it is the right seed
#' for \code{\link{circ_mix}}'s von Mises-family EM. For \eqn{d > 1} the additive
#' distance is a product (independent-coordinate) von Mises model -- a seeding
#' approximation the joint EM density then refines.
#'
#' Starts are chosen by k-means++ on the circular distance (Arthur and
#' Vassilvitskii, 2007), which spreads the initial centres and makes empty
#' clusters rare; the lowest total-within-cluster-distance partition over
#' \code{nstart} starts is returned. The result depends on the random seed, so set
#' one for reproducibility.
#'
#' @param x A numeric vector or matrix of angles in radians, one row per
#'   observation and one column per circular coordinate (a vector is one column).
#' @param centers The number of clusters \eqn{K}.
#' @param nstart The number of k-means++ starts; the lowest-distance partition is
#'   kept.
#' @param iter.max The maximum number of Lloyd iterations per start.
#' @return A list, with the \code{\link[stats]{kmeans}} fields on the circular
#'   metric:
#' \item{cluster}{the length-\code{nrow(x)} integer cluster labels.}
#' \item{centers}{the \eqn{K \times d} matrix of circular-mean centres, in radians
#' on \eqn{(-\pi, \pi]}.}
#' \item{withinss, tot.withinss}{the per-cluster and total within-cluster cosine
#' distance.}
#' \item{size}{the number of points in each cluster.}
#' @references
#' Lloyd, S. P. (1982) Least squares quantization in PCM. \emph{IEEE Transactions
#' on Information Theory} 28, 129-137.
#'
#' Arthur, D. and Vassilvitskii, S. (2007) k-means++: the advantages of careful
#' seeding. \emph{Proceedings of the Eighteenth Annual ACM-SIAM Symposium on
#' Discrete Algorithms}, 1027-1035.
#'
#' Banerjee, A., Dhillon, I. S., Ghosh, J. and Sra, S. (2005) Clustering on the
#' unit hypersphere using von Mises-Fisher distributions. \emph{Journal of Machine
#' Learning Research} 6, 1345-1382.
#'
#' Jammalamadaka, S. R. and SenGupta, A. (2001) \emph{Topics in Circular
#' Statistics}. World Scientific, Singapore.
#' @seealso \code{\link{circ_mix}}, which uses it to seed and split components;
#'   \code{\link[stats]{kmeans}}, the Euclidean original.
#' @examples
#' set.seed(1)
#' theta <- c(stats::rnorm(50, 0, 0.3), stats::rnorm(50, pi, 0.3)) %% (2 * pi)
#' km <- circ_kmeans(theta, 2)
#' km$centers          # two mean directions, near 0 and pi
#' table(km$cluster)
#' @export
circ_kmeans <- function(x, centers, nstart = 5L, iter.max = 50L) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- nrow(x); K <- as.integer(centers)
  if (length(K) != 1L || is.na(K) || K < 1L)
    stop("circ_kmeans: 'centers' must be a single integer >= 1.")
  if (K > n)
    stop("circ_kmeans: more cluster centres (", K, ") than data points (", n, ").")

  one <- function() {
    mu <- x[.circ_kmeans_seed(x, K), , drop = FALSE]   # k-means++ seeds on the torus
    cl <- .circ_assign(x, mu)
    for (it in seq_len(iter.max)) {
      mu  <- .circ_centers(x, cl, K, mu)
      new <- .circ_assign(x, mu)
      if (identical(new, cl)) break
      cl  <- new
    }
    mu <- .circ_centers(x, cl, K, mu)
    D  <- .circ_costs(x, mu)                           # n x K cosine distances
    within <- vapply(seq_len(K),
                     function(k) sum(D[cl == k, k]), numeric(1))
    list(cluster = cl, centers = mu, withinss = within, tot.withinss = sum(within))
  }

  best <- one()
  for (r in seq_len(max(1L, as.integer(nstart)) - 1L)) {
    cand <- one()
    if (cand$tot.withinss < best$tot.withinss) best <- cand
  }
  best$size <- as.integer(tabulate(best$cluster, K))
  best
}

## n x K matrix of summed cosine distances sum_j {1 - cos(x_ij - mu_kj)} -- the
## torus dissimilarity, additive over the d angular coordinates.
.circ_costs <- function(x, mu) {
  out <- matrix(0, nrow(x), nrow(mu))
  for (j in seq_len(ncol(x)))
    out <- out + (1 - cos(outer(x[, j], mu[, j], `-`)))
  out
}

## nearest-centre labels under the circular (cosine) distance.
.circ_assign <- function(x, mu)
  max.col(-.circ_costs(x, mu), ties.method = "first")

## per-cluster circular-mean centres -- the circular mean is the exact minimiser
## of a cluster's summed cosine distance (it maximises sum_i cos(theta_i - mu)), so
## this update is the circular analogue of recentring on the arithmetic mean. An
## empty cluster is reseeded to the currently worst-explained point (largest summed
## cosine distance to its own centre), so a cluster that momentarily empties does
## not collapse for good.
.circ_centers <- function(x, cl, K, mu) {
  for (k in seq_len(K)) {
    ix <- which(cl == k)
    if (length(ix))
      for (j in seq_len(ncol(x)))
        mu[k, j] <- atan2(mean(sin(x[ix, j])), mean(cos(x[ix, j])))
    else
      mu[k, ] <- x[which.max(rowSums(.circ_costs(x, mu))), ]
  }
  mu
}

## k-means++ seed indices (Arthur & Vassilvitskii 2007): the first centre is
## uniform, each next is drawn with probability proportional to D^2, the squared
## Euclidean distance to the nearest centre already chosen. In the (cos, sin)
## embedding D^2 = ||e^{i.theta} - e^{i.mu}||^2 = 2(1 - cos(theta - mu)) = 2 * the
## cosine distance .circ_costs returns -- so the correct D^2 weight IS that cosine
## distance (dmin), NOT its square (squaring would give a D^4 weighting). Well-
## separated starts, far fewer empties.
.circ_kmeans_seed <- function(x, K) {
  n <- nrow(x)
  idx <- integer(K); idx[1L] <- sample.int(n, 1L)
  if (K == 1L) return(idx)
  dmin <- .circ_costs(x, x[idx[1L], , drop = FALSE])[, 1L]   # cosine dist = D^2 / 2
  for (k in 2L:K) {
    idx[k] <- if (sum(dmin) > 0) sample.int(n, 1L, prob = dmin) else sample.int(n, 1L)
    dmin <- pmin(dmin, .circ_costs(x, x[idx[k], , drop = FALSE])[, 1L])
  }
  idx
}

## ---- geometry-aware k-means for circ_mix (init + split) ----------------------
## Cluster a per-unit feature on the response's geometry: circular k-means on the
## circle / torus for an angular response, ordinary stats::kmeans for the linear
## l~c leg. Returns the stats::kmeans-shaped $cluster and $centers either way, so
## the two call sites read the same.
.circ_mix_kmeans <- function(x, K, circular, nstart = 5L) {
  if (circular) circ_kmeans(x, K, nstart = nstart)
  else stats::kmeans(x, centers = K, nstart = nstart)
}

## nearest-centre labels for x given centres, on the response's geometry --
## circular cosine distance for angles, squared Euclidean for the linear leg
## (matching stats::kmeans' own assignment). The split move uses it to label every
## unit by the two sub-centres found on the split component's members.
.circ_mix_assign <- function(x, centers, circular) {
  if (circular) return(.circ_assign(x, centers))
  D <- matrix(0, nrow(x), nrow(centers))
  for (j in seq_len(ncol(x)))
    D <- D + outer(x[, j], centers[, j], `-`)^2
  max.col(-D, ties.method = "first")
}
