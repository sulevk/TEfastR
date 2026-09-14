.OPT_TOL <- 1e-4

.normalize_means <- function(m) {
  s <- sum(m)
  if (s > 0) m / s else m
}

# Vectorized analogue of EMAlgorithm.py::EMUpdate()/computeAbundances(): M is
# a sparse (n_multireads x n_TEs) incidence matrix where M[r, t] counts how
# many of read r's candidate loci overlapped TE t (matching multiplicities in
# the Python flat multi_algn list, including within-read duplicates).
.em_update <- function(means, uniq_counts, eff_len, M) {
  row_mass <- as.numeric(M %*% means)
  inv_row_mass <- ifelse(row_mass > 0, 1 / row_mass, 0)
  col_factor <- as.numeric(Matrix::colSums(Matrix::Diagonal(x = inv_row_mass) %*% M))
  multi_counts <- means * col_factor
  means_out <- ifelse(eff_len > 0, (uniq_counts + multi_counts) / eff_len, 0)
  .normalize_means(means_out)
}

.compute_abundances <- function(means, M) {
  row_mass <- as.numeric(M %*% means)
  inv_row_mass <- ifelse(row_mass > 0, 1 / row_mass, 0)
  means * as.numeric(Matrix::colSums(Matrix::Diagonal(x = inv_row_mass) %*% M))
}

# SQUAREM-accelerated EM (Varadhan & Roland 2008, as used by TEtranscripts'
# EMAlgorithm.py::EMestimate()) for redistributing multi-mapped reads across
# TE instances, proportional to each instance's (uniquely-mapped +
# currently-estimated-multi-mapped) read density.
.em_estimate <- function(uniq_counts, init_multi_counts, M, te_len, num_itr, read_length) {
  eff_len <- te_len - read_length + 1L
  means0 <- .normalize_means(ifelse(eff_len > 0, (uniq_counts + init_multi_counts) / eff_len, 0))

  min_step <- 1; max_step <- 1; m_step <- 4

  cur_iter <- 0L
  while (cur_iter < num_itr) {
    cur_iter <- cur_iter + 1L
    means1 <- .em_update(means0, uniq_counts, eff_len, M)
    means2 <- .em_update(means1, uniq_counts, eff_len, M)

    r  <- means1 - means0
    r2 <- means2 - means1
    v  <- r2 - r

    r_norm <- sqrt(sum(r * r))
    v_norm <- sqrt(sum(v * v))
    r2_norm <- sqrt(sum(r2 * r2))
    rv_norm <- sqrt(abs(sum(r * v)))

    if (v_norm == 0) { means0 <- means1; break }
    alpha <- max(min_step, min(max_step, r_norm / rv_norm))

    if (r_norm < .OPT_TOL) break
    if (r2_norm < .OPT_TOL) { means0 <- means2; break }

    means_prime <- pmax(0, means0 + 2 * alpha * r + (alpha^2) * v)

    if (abs(alpha - 1) > 0.01) {
      means_prime <- .em_update(means_prime, uniq_counts, eff_len, M)
      if (alpha == max_step) { max_step <- max(1, max_step / m_step); alpha <- 1 }
    }
    if (alpha == max_step) max_step <- m_step * max_step

    means0 <- means_prime
  }

  .compute_abundances(means0, M)
}
