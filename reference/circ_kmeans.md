# Circular k-means clustering

k-means on the torus \\T^d\\: a Lloyd iteration whose dissimilarity is
the summed cosine distance \\\sum_j \\1 - \cos(\theta_j - \mu_j)\\\\ and
whose centres are the per-coordinate circular means \\\mu_j =
\mathrm{atan2}(\overline{\sin\theta_j}, \overline{\cos\theta_j})\\. On
angular data this is the right analogue of
[`kmeans`](https://rdrr.io/r/stats/kmeans.html): it respects wrap-around
(\\\theta\\ and \\\theta + 2\pi\\ are the same point), which Euclidean
k-means on the raw radians does not, and its centres stay *on* the
circle (unit resultant) rather than at the radially shrunk arithmetic
mean of the \\(\cos\theta, \sin\theta)\\ embedding, so points are
assigned by true angular distance. For \\d \> 1\\ columns the distance
sums over coordinates, clustering on the product of circles – the same
torus factorisation
[`circ_mix`](https://circstat.github.io/circlss/reference/circ_mix.md)
uses for a joint angular response, which is why `circ_mix` seeds and
splits its components with this routine.

## Usage

``` r
circ_kmeans(x, centers, nstart = 5L, iter.max = 50L)
```

## Arguments

- x:

  A numeric vector or matrix of angles in radians, one row per
  observation and one column per circular coordinate (a vector is one
  column).

- centers:

  The number of clusters \\K\\.

- nstart:

  The number of k-means++ starts; the lowest-distance partition is kept.

- iter.max:

  The maximum number of Lloyd iterations per start.

## Value

A list, with the [`kmeans`](https://rdrr.io/r/stats/kmeans.html) fields
on the circular metric:

- cluster:

  the length-`nrow(x)` integer cluster labels.

- centers:

  the \\K \times d\\ matrix of circular-mean centres, in radians on
  \\(-\pi, \pi\]\\.

- withinss, tot.withinss:

  the per-cluster and total within-cluster cosine distance.

- size:

  the number of points in each cluster.

## Details

Each Lloyd update sets a centre to the circular mean of its members,
which is exactly the minimiser of that cluster's summed cosine distance
– the circular mean direction \\\phi = \mathrm{atan2}(\sum\sin\theta_i,
\sum\cos\theta_i)\\ maximises \\\sum_i \cos(\theta_i - \mu)\\. So the
alternation is coordinate descent on one objective (monotone,
convergent), the circular counterpart of "the arithmetic mean minimises
squared-Euclidean distance" behind ordinary k-means. Equivalently, since
\\1 - \cos(\theta - \mu) = \tfrac12 \\e^{i\theta} - e^{i\mu}\\^2\\, this
is spherical k-means on the unit circle with centres projected back onto
it – the hard-assignment limit of a von Mises mixture with common
concentration (Banerjee et al., 2005), which is why it is the right seed
for
[`circ_mix`](https://circstat.github.io/circlss/reference/circ_mix.md)'s
von Mises-family EM. For \\d \> 1\\ the additive distance is a product
(independent-coordinate) von Mises model – a seeding approximation the
joint EM density then refines.

Starts are chosen by k-means++ on the circular distance (Arthur and
Vassilvitskii, 2007), which spreads the initial centres and makes empty
clusters rare; the lowest total-within-cluster-distance partition over
`nstart` starts is returned. The result depends on the random seed, so
set one for reproducibility.

## References

Lloyd, S. P. (1982) Least squares quantization in PCM. *IEEE
Transactions on Information Theory* 28, 129-137.

Arthur, D. and Vassilvitskii, S. (2007) k-means++: the advantages of
careful seeding. *Proceedings of the Eighteenth Annual ACM-SIAM
Symposium on Discrete Algorithms*, 1027-1035.

Banerjee, A., Dhillon, I. S., Ghosh, J. and Sra, S. (2005) Clustering on
the unit hypersphere using von Mises-Fisher distributions. *Journal of
Machine Learning Research* 6, 1345-1382.

Jammalamadaka, S. R. and SenGupta, A. (2001) *Topics in Circular
Statistics*. World Scientific, Singapore.

## See also

[`circ_mix`](https://circstat.github.io/circlss/reference/circ_mix.md),
which uses it to seed and split components;
[`kmeans`](https://rdrr.io/r/stats/kmeans.html), the Euclidean original.

## Examples

``` r
set.seed(1)
theta <- c(stats::rnorm(50, 0, 0.3), stats::rnorm(50, pi, 0.3)) %% (2 * pi)
km <- circ_kmeans(theta, 2)
km$centers          # two mean directions, near 0 and pi
table(km$cluster)
```
