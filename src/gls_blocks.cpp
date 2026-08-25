// CondPED per-SNP GLS block kernel (Stage: Rcpp migration).
//
// Implementation-layer acceleration ONLY. The statistics are frozen; this
// kernel replicates, chunk-wise and exactly, the R reference
// `.gls_block_components_r()`:
//   U  = Ar %*% x_tilde
//   J1 = sum_j x_tilde_j^2 * A_j        (A_arr:  m x m x n)
//   C  = sum_j x_tilde_j  * AM_j        (AM_arr: m x qm x n)
//   J  = J1 - C %*% G_inv %*% t(C)
// plus the frozen rank / generalised-inverse rule of `.safe_inverse()`
// (symmetric path, relative tolerance rank_tol = sqrt(.Machine$double.eps))
// and the Stage-7.3.5 scale-aware zero floor:
//   max|eigenvalue(J)| <= rank_tol * max|diag(J1)|  =>  rank 0,
//   inverse = zero matrix, condition = Inf, status = "rank_deficient",
//   used_pseudoinverse = TRUE.
// arma::pinv default thresholds are NOT used anywhere.
//
// Status codes returned to R: 0 = "ok", 1 = "rank_deficient", 2 = "failed".

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List gls_blocks_cpp(const arma::mat& X_tilde,   // n x k (rotated dosages)
                          const arma::cube& A_arr,    // m x m x n (Vinv)
                          const arma::cube& AM_arr,   // m x qm x n (A_j M_j)
                          const arma::mat& Ar,        // m x n (A_j r_j)
                          const arma::mat& G_inv,     // qm x qm (XtVinvX)^{-1}
                          double rank_tol) {
  const arma::uword n = X_tilde.n_rows;
  const arma::uword k = X_tilde.n_cols;
  const arma::uword m = Ar.n_rows;
  const arma::uword qm = G_inv.n_rows;
  const double NA_D = NA_REAL;

  arma::mat U_mat(m, k), beta_mat(m, k);
  arma::cube J_cube(m, m, k), Jinv_cube(m, m, k);
  arma::ivec rank(k), status(k);
  arma::vec condition(k);
  arma::uvec used_pseudo(k);

  for (arma::uword s = 0; s < k; ++s) {
    const arma::vec x = X_tilde.col(s);
    const arma::vec U = Ar * x;

    arma::mat J1(m, m, arma::fill::zeros);
    arma::mat C(m, qm, arma::fill::zeros);
    {
      // pointer loops over contiguous cube slices: no temporaries, and the
      // element order (j outer, column-major inner) matches the R reference
      // accumulation exactly
      double* j1p = J1.memptr();
      double* cp = C.memptr();
      const size_t mm = m * m, mqm = m * qm;
      for (arma::uword j = 0; j < n; ++j) {
        const double xj = x[j];
        const double x2 = xj * xj;
        const double* aj = A_arr.slice_memptr(j);
        const double* amj = AM_arr.slice_memptr(j);
        for (size_t t = 0; t < mm; ++t) j1p[t] += x2 * aj[t];
        for (size_t t = 0; t < mqm; ++t) cp[t] += xj * amj[t];
      }
    }
    arma::mat J = J1 - C * G_inv * C.t();

    // ---- frozen rank / generalised-inverse rule (eigen path) ----------
    arma::vec ev;
    arma::mat evec;
    arma::mat J_inv(m, m);
    arma::uword r = 0;
    double cond = NA_D;
    int st = 2;              // failed
    unsigned int pseudo = 0; // NA by default (failure)

    if (arma::eig_sym(ev, evec, J)) {
      const double max_abs = arma::abs(ev).max();
      const arma::uvec keep = arma::find(arma::abs(ev) > rank_tol * max_abs);

      if (max_abs == 0.0) {
        // zero spectrum: MP inverse of a zero matrix is the zero matrix
        J_inv.zeros();
        r = 0;
        cond = arma::datum::inf;
        st = 1;              // rank_deficient
        pseudo = 1;
      } else {
        const double cutoff = rank_tol * max_abs;
        r = keep.n_elem;
        const double min_ev = ev.min();   // eig_sym: ascending order
        if (r == m && min_ev > cutoff) {
          // numerically positive definite: Cholesky fast path
          // (J = R' R with upper R; J^{-1} = R^{-1} (R^{-1})', the
          // chol2inv computation)
          arma::mat Rc;
          if (arma::chol(Rc, J)) {
            const arma::mat Ri = arma::inv(arma::trimatu(Rc));
            J_inv = Ri * Ri.t();
            cond = ev.max() / min_ev;
            st = 0;          // ok
            pseudo = 0;
          } else {
            // fall through to the eigen pseudo-inverse (full rank, ok)
            const arma::vec d_inv = 1.0 / ev;
            J_inv = evec * arma::diagmat(d_inv) * evec.t();
            cond = max_abs / arma::abs(ev).min();
            st = 0;
            pseudo = 1;
          }
        } else {
          arma::vec d_inv(m, arma::fill::zeros);
          d_inv.elem(keep) = 1.0 / ev.elem(keep);
          J_inv = evec * arma::diagmat(d_inv) * evec.t();
          cond = (r < m) ? arma::datum::inf
                         : max_abs / arma::abs(ev).min();
          st = (r < m) ? 1 : 0;
          pseudo = 1;
        }

        // ---- scale-aware zero floor (Stage 7.3.5) ----------------------
        const double J1_scale = arma::abs(J1.diag()).max();
        if (J1_scale > 0.0 && max_abs <= rank_tol * J1_scale) {
          J_inv.zeros();
          r = 0;
          cond = arma::datum::inf;
          st = 1;            // rank_deficient
          pseudo = 1;
        }
      }
    } else {
      J_inv.fill(NA_D);
      r = 0;
      cond = arma::datum::inf;
      st = 2;                // failed
      pseudo = 0;
    }

    U_mat.col(s) = U;
    J_cube.slice(s) = J;
    Jinv_cube.slice(s) = J_inv;
    rank[s] = static_cast<int>(r);
    status[s] = st;
    condition[s] = cond;
    used_pseudo[s] = pseudo;
    beta_mat.col(s) = (r > 0 && st != 2) ? J_inv * U
                                         : arma::vec(m).fill(NA_D);
  }

  return Rcpp::List::create(
    Rcpp::Named("U") = U_mat,
    Rcpp::Named("J") = J_cube,
    Rcpp::Named("J_inv") = Jinv_cube,
    Rcpp::Named("rank") = rank,
    Rcpp::Named("condition") = condition,
    Rcpp::Named("status") = status,
    Rcpp::Named("used_pseudoinverse") = used_pseudo,
    Rcpp::Named("beta") = beta_mat
  );
}
