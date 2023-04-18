#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



// [[Rcpp::export]]

Rcpp::List location_search(int num_times, int num_subjects,
                           arma::cube W_total) {
  int count = 0;
  for (int i = 0; i < num_times; i++) {
    for (int j = 0; j < num_times; j++) {
      for (int s = 0; s < num_subjects; s++) {
        if (W_total(i, j, s) != 0){
          count++;
        }
      }
    }
  }
  int idx = 0;
  arma::vec y(count);
  arma::vec row_location(count);
  arma::vec column_location(count);
  for (int i = 0; i < num_times; i++) {
    for (int j = 0; j < num_times; j++) {
      for (int s = 0; s < num_subjects; s++) {
        if (W_total(i, j, s) != 0) {
          y(idx) = W_total(i, j, s);
          row_location(idx) = i;
          column_location(idx) = j;
          idx++;
        }
      }
    }
  }

  Rcpp::List location;
  location["y"] = y;
  location["row_location"] = row_location;
  location["column_location"] = column_location;

  return location;
}



// [[Rcpp::export]]

arma::cube V_projection(int K, int num_times, int num_subjects,
                        arma::cube X, arma::mat V0) {
  arma::mat Yt(num_subjects, num_times);
  arma::cube Y(num_subjects, num_times, K);
  arma::mat Ytmp(num_subjects, K);
  for (int k = 0; k < K; k++) {
    for (int t = 0; t < num_times; t++) {
      Ytmp = X.col(t);
      Yt.col(t) = Ytmp * V0.col(k);
    }
    Y.slice(k) = Yt;
  }
  return(Y);
}



// [[Rcpp::export]]

arma::cube empirical_W_total(int num_times, int num_subjects, int K,
                             arma::cube OBS, arma::cube Y) {
  arma::cube W_total(num_times, num_times, num_subjects);
  W_total.zeros(num_times, num_times, num_subjects);
  arma::mat OBStmp(num_subjects, num_times);
  arma::mat Ytmp(num_times, K);
  OBStmp = OBS.slice(0);
  for (int s = 0; s < num_times; s++) {
    for (int t = 0; t < num_times; t++) {
      if (s <= t){
        for (int i = 0; i < num_subjects; i++) {
          if (OBStmp(i, t) * OBStmp(i, s) == 1) {
            Ytmp = Y.row(i);
            W_total(s, t, i)= sum(Ytmp.row(s) % Ytmp.row(t));
          }
        }
      }
    }
  }
  for (int i = 0; i < num_times; i++) {
    for (int t = 0; t < num_times; t++) {
      for (int j = 0; j < num_subjects; j++) {
        if (fabs(W_total(i, t, j)) < 1e-10) {
          W_total(i, t, j) = 0.0;
        }
      }
    }
  }
  return(W_total);
}



// [[Rcpp::export]]

int pairs_count(int num_times, int num_subjects, arma::cube W_total) {
  int count = 0;
  arma::mat W_count(num_times, num_times);
  W_count.zeros(num_times, num_times);
  arma::mat W_totaltmp(num_times, num_times);
  for (int i = 0; i < num_subjects; i++) {
    W_totaltmp = W_total.slice(i);
    W_count += W_totaltmp;
  }
  for (int s = 0; s < num_times; s++) {
    for (int t = 0; t < num_times; t++) {
      if (fabs(W_count(s, t)) > 0.0) {
        count++;
      }
    }
  }
  return(count);
}



// [[Rcpp::export]]

arma::mat ridge_penalized(int num_subjects, int num_times, int num_features,
                          int K, arma::mat XI, arma::mat OI, arma::mat V0Phi0,
                          double eps) {
  arma::vec count(num_subjects);
  count.zeros(num_subjects);
  arma::vec X0(num_times * num_features);
  arma::vec O0(num_times * num_features);
  for (int i = 0; i < num_subjects; i++) {
    O0 = trans(OI.row(i));
    for (int j = 0; j < (num_times * num_features); j++) {
      if (O0(j) == 1) {
        count(i)++;
      }
    }
  }

  int c = 0;
  arma::mat a(K * K, K * K);
  arma::vec b(K * K);
  arma::mat Utmp(num_subjects, K * K);
  arma::mat Kdiag(K * K, K * K);
  for (int i = 0; i < num_subjects; i++) {
    X0 = trans(XI.row(i));
    O0 = trans(OI.row(i));
    arma::vec Oidx(count(i));
    arma::vec X01(count(i));
    arma::mat V0Phi01(count(i), K * K);
    int index = 0;
    for (int j = 0; j < (num_times * num_features); j++) {
      if (O0(j) == 1) {
        Oidx(index) = j;
        index++;
      }
    }
    for (int s = 0; s < count(i); s++) {
      c = Oidx(s);
      V0Phi01.row(s) = V0Phi0.row(c);
      X01(s) = X0(c);
    }
    a = trans(V0Phi01) * V0Phi01 + eps * Kdiag.eye(K * K, K * K);
    b = trans(V0Phi01) * X01;
    Utmp.row(i) = trans(solve(a, b));
  }

  return (Utmp);
}



// [[Rcpp::export]]

arma::cube rearrange_G(int K, arma::mat Utmp, arma::mat U0, arma::cube G) {
  arma::mat tmp2(K, K);
  tmp2 = trans(U0) * Utmp;
  arma::vec index(K);
  for (int i = 0; i < K; i++) {
    for (int j = 0; j < K; j++) {
      index(0) = j * K;
      for (int s = 1; s < K; s++) {
        index(s) = index(s - 1) + 1;
      }
      for (int t = 0; t < K; t++) {
        G(i, t, j) = tmp2(i, index(t));
      }
    }
  }
  return(G);
}



// [[Rcpp::export]]

Rcpp::List kronecker_VPhi(int K, arma::mat V0, arma::mat Phi0,
                          arma::mat V, arma::mat Phi, arma::mat VPhi,
                          arma::mat U2, arma::mat U3) {
  for (int i = 0; i < K; i++){
    V.col(i) = V0 * U3.col(i);
    Phi.col(i) = Phi0 * U2.col(i);
    VPhi.col(i) = trans(arma::kron(trans(V.col(i)), trans(Phi.col(i))));
  }

  Rcpp::List VPhi_list;
  VPhi_list["V"] = V;
  VPhi_list["Phi"] = Phi;
  VPhi_list["VPhi"] = VPhi;

  return(VPhi_list);
}


// [[Rcpp::export]]

Rcpp::List rescale_VPhi(int K, double Phi_norm, double V_norm,
                        arma::mat V, arma::mat Phi, arma::mat VPhi) {
  double c2;
  double c3;
  for (int i = 0; i < K; i++) {
    c2 = sqrt(sum(square(Phi.col(i)))) / Phi_norm;
    c3 = sqrt(sum(square(V.col(i)))) / V_norm;
    Phi.col(i) = Phi.col(i) / c2;
    V.col(i) = V.col(i) / c3;
    VPhi.col(i) = trans(arma::kron(trans(V.col(i)), trans(Phi.col(i))));
  }

  Rcpp::List rescale_VPhi_list;
  rescale_VPhi_list["V"] = V;
  rescale_VPhi_list["Phi"] = Phi;
  rescale_VPhi_list["VPhi"] = VPhi;

  return(rescale_VPhi_list);
}



// [[Rcpp::export]]

Rcpp::List rescale_U(int num_subjects, int num_times, int num_features, int K,
                     arma::mat XI, arma::mat OI, arma::mat VPhi, double eps) {
  arma::vec count(num_subjects);
  count.zeros(num_subjects);
  arma::vec x(num_times * num_features);
  arma::vec o(num_times * num_features);
  for (int i = 0; i < num_subjects; i++) {
    o = trans(OI.row(i));
    for (int j = 0; j < (num_times * num_features); j++) {
      if (o(j) == 1) {
        count(i)++;
      }
    }
  }

  int c = 0;
  int d = 0;
  arma::mat a(K, K);
  arma::vec b(K);
  arma::mat U(num_subjects, K);
  arma::mat Kdiag(K, K);
  arma::mat xhat(num_subjects, num_times * num_features);
  xhat.zeros(num_subjects, num_times * num_features);
  for (int i = 0; i < num_subjects; i++) {
    x = trans(XI.row(i));
    o = trans(OI.row(i));
    arma::vec Oidx(count(i));
    arma::vec x1(count(i));
    arma::mat VPhi1(count(i), K);
    int index = 0;
    for (int j = 0; j < (num_times * num_features); j++){
      if (o(j) == 1) {
        Oidx(index) = j;
        index++;
      }
    }
    for (int s = 0; s < count(i); s++) {
      c = Oidx(s);
      VPhi1.row(s) = VPhi.row(c);
      x1(s) = x(c);
    }
    a = trans(VPhi1) * VPhi1 + eps * Kdiag.eye(K, K);
    b = trans(VPhi1) * x1;
    U.row(i) = trans(solve(a, b));
    for (int t = 0; t < count(i); t++) {
      d = Oidx(t);
      xhat(i, d) = sum(VPhi.row(d) % U.row(i));
    }
  }

  Rcpp::List rescale_U_list;
  rescale_U_list["U"] = U;
  rescale_U_list["xhat"] = xhat;
  return (rescale_U_list);
}



// [[Rcpp::export]]

arma::vec sigma_init(bool homoNoise, int num_features, int Osum,
                     arma::mat xhat, arma::mat XI, arma::mat index) {
  arma::vec sigma(num_features);
  arma::vec xhattmp(Osum);
  arma::vec XItmp(Osum);
  if (homoNoise){
    for (int j = 0; j < num_features; j++) {
      for (int i = 0; i < Osum; i++) {
        xhattmp(i) = xhat(index(i,0), index(i,1));
        XItmp(i) = XI(index(i,0), index(i,1));
      }
      sigma(j) = sum(square(xhattmp - XItmp)) / Osum;
    }
  }
  return(sigma);
}



// [[Rcpp::export]]

arma::mat posterior_kronecker_VPhi(int K, arma::mat V, arma::mat Phi,
                                   arma::mat VPhi) {
  for (int i = 0; i < K; i++) {
    VPhi.col(i) = arma::kron(V.col(i), Phi.col(i));
  }
  return(VPhi);
}



// [[Rcpp::export]]

arma::mat sigma_transform(int num_times, int num_features,
                          arma::vec sigma, arma::vec s) {
  arma::vec idx(num_times);
  for (int i = 0; i < num_features; i++) {
    idx(0) = i * num_times;
    for (int j = 1; j < num_times; j++) {
      idx(j) = idx(j - 1) + 1;
    }
    for (int t = 0; t < num_times; t++) {
      s(idx(t)) = sigma(i);
    }
  }
  return (s);
}



// [[Rcpp::export]]

Rcpp::List posterior_distribution(int num_subjects, int num_times,
                                  int num_features, int K, bool fit_intercept,
                                  arma::vec intercepts, arma::mat XI,
                                  arma::mat OI, arma::mat VPhi, arma::vec s,
                                  arma::vec sigmaF, arma::mat Z,
                                  arma::mat beta) {
  arma::vec count(num_subjects);
  count.zeros(num_subjects);
  arma::vec x1(num_times * num_features);
  arma::vec o1(num_times * num_features);
  for (int i = 0; i < num_subjects; i++) {
    o1 = trans(OI.row(i));
    for (int j = 0; j < (num_times * num_features); j++) {
      if (o1(j) == 1) {
        count(i)++;
      }
    }
  }

  int c = 0;
  arma::mat vec(num_subjects, K);
  vec.zeros(num_subjects, K);
  arma::cube mat(num_subjects, K, K);
  mat.zeros(num_subjects, K, K);
  arma::mat mu(num_subjects, K);
  mu.zeros(num_subjects, K);
  arma::cube cov(num_subjects, K, K);
  cov.zeros(num_subjects, K, K);
  for (int i = 0; i < num_subjects; i++) {
    x1 = trans(XI.row(i));
    o1 = trans(OI.row(i));
    arma::vec Oidx(count(i));
    arma::vec x10(count(i));
    arma::vec s0(count(i));
    arma::mat VPhi0(count(i), K);
    int index = 0;
    for (int j = 0; j < (num_times * num_features); j++) {
      if (o1(j) == 1) {
        Oidx(index) = j;
        index++;
      }
    }
    for (int t = 0; t < count(i); t++) {
      c = Oidx(t);
      VPhi0.row(t) = VPhi.row(c);
      x10(t) = x1(c);
      s0(t) = s(c);
    }

    arma::mat tmp1(K, K);
    arma::vec tmp2(K);
    arma::vec tmp3(K);
    arma::mat tmp4(K, K);
    arma::mat M1(K, K);
    arma::mat M2(K, K);
    arma::mat covtmp(K, K);
    arma::mat VPhi0tmp(K, count(i));
    arma::mat VPhi0trans(K, count(i));
    M1.zeros(K, K);
    VPhi0trans = trans(VPhi0);
    for (int n = 0; n < K; n++) {
      VPhi0tmp.row(n) = VPhi0trans.row(n) / trans(s0);
    }
    tmp1 = VPhi0tmp * VPhi0;
    tmp2 = trans(VPhi0) * (x10 / s0);
    for (int m = 0; m < K; m++) {
      M1(m, m) = 1.0 / sigmaF(m);
    }
    tmp4 = tmp1 + M1;
    cov.row(i) =  inv(tmp4);
    covtmp = cov.row(i);
    M2 = M1 - (M1 * covtmp) * M1;
    tmp3.zeros(K);
    tmp3 = trans(Z.row(i) * beta) / sigmaF;
    if (fit_intercept) {
      tmp3 += intercepts / sigmaF;
    }
    mu.row(i) = trans(covtmp * (tmp2 + tmp3));
    mat.row(i) = M2;
    vec.row(i) = (trans(tmp2) * covtmp) * M1;
  }

  Rcpp::List posterior;
  posterior["mu"] = mu;
  posterior["cov"] = cov;
  posterior["mat"] = mat;
  posterior["vec"] = vec;

  return(posterior);
}



// [[Rcpp::export]]

arma::mat kronecker_Vmu(int K, int num_subjects, int num_features,
                        arma::mat V, arma::mat mu, arma::mat A0) {
  for (int i = 0; i < K; i++) {
    A0.col(i) = arma::kron(V.col(i), mu.col(i));
  }
  return(A0);
}



// [[Rcpp::export]]

arma::vec s2_transform(int num_subjects, int num_features,
                       arma::vec sigma, arma::vec s2) {
  arma::vec idx(num_features);
  for (int i = 0; i < num_subjects; i++) {
    idx(0) = i * num_features;
    for (int j = 1; j < num_features; j++) {
      idx(j) = idx(j - 1) + 1;
    }
    for (int t = 0; t < num_features; t++) {
      s2(idx(t)) = sigma(i);
    }
  }
  return(s2);
}



// [[Rcpp::export]]

Rcpp::List Mt_creator(int num_subjects, int num_times, int num_features, int K,
                      arma::vec s2, arma::vec sigma, arma::mat XT, arma::mat OT,
                      arma::mat A0, arma::mat V,
                      arma::cube cov, arma::cube OBS) {
  arma::vec x0(num_subjects * num_features);
  arma::vec A1(num_subjects * num_features);
  arma::mat h(K, num_times);
  arma::vec x00(num_subjects * num_features);
  arma::vec A11(num_subjects * num_features);
  arma::vec A2(num_subjects);
  arma::vec A3(num_features);
  arma::vec A4(num_features);
  arma::mat covtmp(num_subjects, K);
  arma::mat OBStmp1(num_subjects, num_times);
  arma::vec OBStmp2(num_subjects);
  arma::cube M(K, K, num_times);
  for (int k = 0; k < K; k++) {
    for (int t = 0; t < num_times; t++) {
      x0 = trans(XT.row(t));
      A1 = A0.col(k);
      int count1 = 0;
      arma::vec OTtmp(num_subjects * num_features);
      OTtmp = trans(OT.row(t));
      for (int i = 0; i < (num_subjects * num_features); i++) {
        if (OTtmp(i) == 1) {
          count1++;
        }
      }
      arma::vec idx1(count1);
      int index1 = 0;
      for (int p = 0; p < (num_subjects * num_features); p++) {
        if (OTtmp(p) == 1) {
          idx1(index1) = p;
          index1++;
        }
      }
      arma::vec x0tmp(count1);
      arma::vec A1tmp(count1);
      arma::vec s2tmp(count1);
      for (int q = 0; q < count1; q++) {
        x0tmp(q) = x0(idx1(q));
        A1tmp(q) = A1(idx1(q));
        s2tmp(q) = s2(idx1(q));
      }
      h(k, t) = sum(x0tmp % A1tmp / s2tmp);
      for (int s = 0; s < K; s++) {
        if (s <= k) {
          x00 = A0.col(k);
          A11 = A0.col(s);
          covtmp = cov.slice(s);
          A2 = covtmp.col(k);

          arma::vec x00tmp(count1);
          arma::vec A11tmp(count1);
          arma::vec s22tmp(count1);
          for (int q = 0; q < count1; q++) {
            x00tmp(q) = x00(idx1(q));
            A11tmp(q) = A11(idx1(q));
            s22tmp(q) = s2(idx1(q));
          }
          M(k, s, t) = sum(x00tmp % A11tmp / s22tmp);

          for (int j = 0; j < num_features; j++) {
            OBStmp1 = OBS.slice(j);
            OBStmp2 = OBStmp1.col(t);
            int count2 = 0;
            for (int n = 0; n < num_subjects; n++) {
              if (OBStmp2(n) == 1) {
                count2++;
              }
            }
            arma::vec idx2(count2);
            int index2 = 0;
            for (int l = 0; l < num_subjects; l++) {
              if (OBStmp2(l) == 1) {
                idx2(index2) = l;
                index2++;
              }
            }
            arma::vec array(num_subjects);
            for (int m = 0; m < num_subjects; m++) {
              array(m) = m;
            }
            arma::vec idx3(count2);
            int index3 = 0;
            for (int z = 0; z < count2; z++) {
              index3 = idx2(z);
              idx3(z) = array(index3);
            }

            arma::vec A2tmp(count2);
            arma::vec s222tmp(count2);
            for (int q = 0; q < count2; q++) {
              A2tmp(q) = A2(idx2(q));
              s222tmp(q) = s2(idx2(q));
            }
            double V1 = V(j, k);
            double V2 = V(j, s);
            M(k, s, t) += sum(A2tmp * V1 * V2 / s222tmp);
          }
          M(s, k, t) = M(k, s, t);
        }
      }
    }
  }

  Rcpp::List Mt_list;
  Mt_list["M"] = M;
  Mt_list["h"] = h;

  return(Mt_list);
}



// [[Rcpp::export]]

arma::vec penalty_search(int num_times, int nlambda1, int max_iter0,
                         arma::vec lambda1s, arma::mat A, arma::mat Omega,
                         arma::vec dfs) {
  int iter0;
  double tol0;
  double error0;
  arma::vec tmp1(num_times);
  arma::vec tmp2(num_times);
  arma::mat U1(num_times, num_times);
  arma::mat V1(num_times, num_times);
  arma::mat U2(num_times, num_times);
  arma::mat V2(num_times, num_times);
  double maxtmp;
  double mintmp;
  double df0;
  double middletmp = 0.0;
  double dfhat;
  arma::mat invA(num_times, num_times);
  arma::vec diagA(num_times);
  arma::mat diagAtmp(num_times, num_times);
  for (int j = 0; j < nlambda1; j++) {
    iter0 = 0;
    tol0 = 1e-4;
    error0 = tol0 * 2;
    svd(U1, tmp1, V1, A);
    svd(U2, tmp2, V2, Omega);
    maxtmp = max(tmp1) / (min(tmp2) + 1e-6);
    mintmp = 1e-10;
    df0 = dfs(j);
    while (iter0 < max_iter0 && fabs(error0) > tol0) {
      middletmp = (maxtmp + mintmp) / 2;
      invA = inv(A + middletmp * Omega);
      diagAtmp = invA * A;
      for (int t = 0; t < num_times; t++) {
        diagA(t) = diagAtmp(t, t);
      }
      dfhat = sum(diagA);
      error0 = dfhat - df0;
      if (error0 > tol0) {
        mintmp = middletmp;
      }
      else if (error0 < -tol0) {
        maxtmp = middletmp;
      }
      iter0++;
    }
    lambda1s(j) = middletmp;
  }
  return(lambda1s);
}



// [[Rcpp::export]]

arma::mat Omega_errors(int num_times, int nlambda1, double ridge_traj,
                       arma::vec a, arma::vec b, arma::vec lambda1s,
                       arma::vec phi, arma::mat Omega, arma::mat errors) {
  arma::mat tmp(num_times, num_times);
  arma::vec ridge(num_times);
  ridge.fill(ridge_traj);
  arma::vec a_ridge(num_times);
  a_ridge = a + ridge;
  arma::mat a_ridge_diag(num_times, num_times);
  a_ridge_diag.eye(num_times, num_times);
  for (int l = 0; l < num_times; l++) {
    a_ridge_diag(l, l) = a_ridge(l);
  }
  arma::vec phi0(num_times);
  arma::vec H_diag(num_times);
  for (int s = 0; s < nlambda1; s++) {
    tmp = a_ridge_diag + lambda1s(s) * Omega;
    phi0 = solve(tmp, b);
    tmp = inv(tmp);
    for (int t = 0; t < num_times; t++) {
      H_diag(t) = tmp(t, t);
    }
    H_diag %= a_ridge;
    errors.col(s) = a_ridge % square(phi0 - phi) / square(1.0 - H_diag);
  }
  return errors;
}



// [[Rcpp::export]]

double Omega_update(int num_times, int nlambda1, int lambda1_dfmin,
                    int lambda1_dfmax, double ridge_traj, double lam1criterion,
                    arma::vec a, arma::vec b, arma::mat Omega) {
  arma::mat A(num_times, num_times);
  for (int t = 0; t < num_times; t++) {
    A(t, t) = a(t);
  }
  arma::vec dfs_tmp(nlambda1);
  for (int l = 0; l < nlambda1; l++) {
    dfs_tmp(l) = 1.0 - l / (nlambda1 - 1.0);
  }
  arma::vec dfs(nlambda1);
  dfs = lambda1_dfmin + dfs_tmp * (lambda1_dfmax - 1.0);
  arma::vec lambda1s(num_times);
  lambda1s.zeros(num_times);
  lambda1s = penalty_search(num_times, nlambda1, 100, lambda1s, A, Omega, dfs);
  arma::vec ridge(num_times);
  ridge.fill(ridge_traj);
  arma::vec a_ridge(num_times);
  a_ridge = a + ridge;
  arma::vec phi(num_times);
  phi = b / a_ridge;

  arma::mat errors(num_times, nlambda1);
  arma::mat tmp(num_times, num_times);
  arma::mat a_ridge_diag(num_times, num_times);
  a_ridge_diag.eye(num_times, num_times);
  for (int l = 0; l < num_times; l++) {
    a_ridge_diag(l, l) = a_ridge(l);
  }
  arma::vec phi0(num_times);
  arma::vec H_diag(num_times);
  for (int s = 0; s < nlambda1; s++) {
    tmp = a_ridge_diag + lambda1s(s) * Omega;
    phi0 = solve(tmp, b);
    tmp = inv(tmp);
    for (int t = 0; t < num_times; t++) {
      H_diag(t) = tmp(t, t);
    }
    H_diag %= a_ridge;
    errors.col(s) = a_ridge % square(phi0 - phi) / square(1.0 - H_diag);
  }

  arma::vec errors_mean(nlambda1);
  for (int j = 0; j < nlambda1; j++) {
    double col_sum = 0.0;
    for (int i = 0; i < num_times; i++) {
      col_sum += errors(i, j);
    }
    errors_mean(j) = col_sum / num_times;
  }

  arma::vec errors_sd(nlambda1);
  for (int j = 0; j < nlambda1; j++) {
    double col_sum_squares = 0.0;
    for (int i = 0; i < num_times; i++) {
      col_sum_squares += (errors(i, j) - errors_mean(j)) *
        (errors(i, j) - errors_mean(j)) / num_times;
    }
    errors_sd(j) = sqrt(col_sum_squares / num_times);
  }

  int idx1;
  idx1 = index_min(errors_mean);
  double c;
  c = errors_mean(idx1) + errors_sd(idx1) * lam1criterion;
  int count = 0;
  for (int m = 0; m < nlambda1; m++) {
    if (errors_mean(m) <= c) {
      count ++;
    }
  }
  arma::vec idx2s(count);
  int idx = 0;
  for (int n = 0; n < nlambda1; n++) {
    if (errors_mean(n) <= c) {
      idx2s(idx) = n;
      idx++;
    }
  }
  int idx2;
  idx2 = max(idx2s);
  double lambda1;
  lambda1 = lambda1s(idx2);
  return lambda1;
}



// [[Rcpp::export]]

double binary_search(int num, double norm_constraint,
                     arma::vec a, arma::vec b) {
  int max_iter = 1000;
  if (min(a) <= 0) {
    Rcpp::stop("Encountered nonPD input during the subroutine updating V.");
  }
  arma::vec z(num);
  z = b / a;
  double c;
  c = sum(b % b);
  double norm2_0;
  norm2_0 = sqrt(sum(square(z)));
  double lam_min;
  double lam_max;
  int idx;
  idx = index_min(a);
  if (norm2_0 > norm_constraint + 1e-8) {
    lam_min = 0.0;
    lam_max = sqrt(c) / norm_constraint - min(a);
  } else if (norm2_0 < norm_constraint - 1e-8) {
    lam_min = -(a(idx) - fabs(b(idx)) / norm_constraint);
    lam_max = 0.0;
  } else {
    lam_min = 0.0;
    lam_max = 0.0;
  }

  double lam_current = 0.0;
  int iter = 0;
  while (iter < max_iter && fabs(norm2_0 - norm_constraint) > 1e-8) {
    arma::vec current(num);
    if (norm2_0 > norm_constraint) {
      lam_min = lam_current;
      lam_current = (lam_current + lam_max) / 2.0;
    } else {
      lam_max = lam_current;
      lam_current = (lam_current + lam_min) / 2.0;
    }
    current.fill(lam_current);
    z = b / (a + current);
    norm2_0 = sqrt(sum(square(z)));
    iter++;
  }

  return lam_current;
}



// [[Rcpp::export]]

arma::vec phi_transform(int num_times, double lambda1, double ridge_traj,
                        double h, arma::vec a, arma::vec b, arma::mat Omega) {
  arma::vec ridge(num_times);
  ridge.fill(ridge_traj);
  arma::vec a_ridge(num_times);
  a_ridge = a + ridge;
  arma::mat a_ridge_diag(num_times, num_times);
  a_ridge_diag.eye(num_times, num_times);
  for (int l = 0; l < num_times; l++) {
    a_ridge_diag(l, l) = a_ridge(l);
  }
  arma::mat A(num_times, num_times);
  A = Omega * lambda1 + a_ridge_diag;
  arma::vec tmp(num_times);
  arma::mat U(num_times, num_times);
  arma::mat V(num_times, num_times);
  svd(U, tmp, V, A);
  arma::vec a1(num_times);
  arma::vec b1(num_times);
  a1 = tmp;
  b1 = trans(V) * b;
  double laplace_lam;
  laplace_lam = binary_search(num_times, h, a1, b1);
  arma::vec phi(num_times);
  arma::vec lam(num_times);
  lam.fill(laplace_lam);
  phi = b1 / (a1 + lam);
  phi = U * phi;
  return phi;
}



// [[Rcpp::export]]

Rcpp::List Phi_solver(int K, int num_times, int nlambda1,
                      int lambda1_dfmin, int lambda1_dfmax,
                      double ridge_traj, double lam1criterion,
                      double self_h, bool update_smooth_penalty,
                      arma::vec lambda1, arma::mat self_Phi, arma::mat h,
                      arma::mat Omega, arma::cube M) {
  for (int k = 0; k < K; k++) {
    arma::mat Mk(K, num_times);
    arma::vec a(num_times);
    arma::vec b(num_times);
    Mk = M.row(k);
    a = trans(Mk.row(k));
    b.zeros(num_times);
    b += trans(h.row(k));
    for (int s = 0; s < K; s++) {
      if (s != k) {
        b -= trans(Mk.row(s)) % self_Phi.col(s);
      }
    }
    arma::vec phi(num_times);
    if (update_smooth_penalty) {
      lambda1(k) = Omega_update(num_times, nlambda1, lambda1_dfmin,
              lambda1_dfmax, ridge_traj, lam1criterion, a, b, Omega);
      phi = phi_transform(num_times, lambda1(k), ridge_traj, sqrt(self_h),
                            a, b, Omega);
    } else {
      phi = phi_transform(num_times, 0.0, ridge_traj, sqrt(self_h),
                            a, b, Omega);
    }
    self_Phi.col(k) = phi;
  }

  Rcpp::List Phi_solver_list;
  Phi_solver_list["Phi"] = self_Phi;
  Phi_solver_list["lambda1"] = lambda1;
  return Phi_solver_list;
}



// [[Rcpp::export]]

arma::mat kronecker_Phimu(int K, int num_subjects, int num_times,
                          arma::mat Phi, arma::mat mu, arma::mat A0){
  for (int i = 0; i < K; i++){
    A0.col(i) = arma::kron(Phi.col(i), mu.col(i));
  }
  return(A0);
}



// [[Rcpp::export]]

Rcpp::List Mj_creator(int num_subjects, int num_times, int num_features, int K,
                      arma::vec sigma, arma::mat XJ, arma::mat OJ, arma::mat A0,
                      arma::mat Phi, arma::cube cov, arma::cube OBS){
  arma::vec x0(num_subjects * num_times);
  arma::vec A1(num_subjects * num_times);
  arma::mat h(K, num_features);
  arma::vec A2(num_subjects);
  arma::mat covtmp(num_subjects, K);
  arma::cube M(K, K, num_features);
  arma::mat OBStmp1(num_subjects, num_times);
  arma::vec OBStmp2(num_subjects);
  for (int k = 0; k < K; k++){
    for (int j = 0; j < num_features; j++){
      x0 = trans(XJ.row(j));
      A1 = A0.col(k);
      int count1 = 0;
      arma::vec OJtmp(num_subjects * num_times);
      OJtmp = trans(OJ.row(j));
      for (int i = 0; i < (num_subjects * num_times); i++){
        if (OJtmp(i) == 1){
          count1++;
        }
      }
      arma::vec idx1(count1);
      int index1 = 0;
      for (int p = 0; p < (num_subjects * num_times); p++){
        if (OJtmp(p) == 1){
          idx1(index1) = p;
          index1++;
        }
      }
      arma::vec x0tmp(count1);
      arma::vec A1tmp(count1);
      for (int q = 0; q < count1; q++){
        x0tmp(q) = x0(idx1(q));
        A1tmp(q) = A1(idx1(q));
      }
      h(k, j) = sum(x0tmp % A1tmp) / sigma(j);
      for (int s = 0; s < K; s++){
        if (s <= k){
          x0 = A0.col(k);
          A1 = A0.col(s);
          covtmp = cov.slice(s);
          A2 = covtmp.col(k);
          for (int m = 0; m < count1; m++){
            x0tmp(m) = x0(idx1(m));
            A1tmp(m) = A1(idx1(m));
          }
          M(k, s, j) = sum(x0tmp % A1tmp);
          for (int t = 0; t < num_times; t++){
            OBStmp1 = OBS.slice(j);
            OBStmp2 = OBStmp1.col(t);
            int count2 = 0;
            for (int n = 0; n < num_subjects; n++){
              if (OBStmp2(n) == 1){
                count2++;
              }
            }
            arma::vec idx2(count2);
            int index2 = 0;
            for (int l = 0; l < num_subjects; l++){
              if (OBStmp2(l) == 1){
                idx2(index2) = l;
                index2++;
              }
            }
            arma::vec idx3(num_subjects);
            for (int w = 0; w < num_subjects; w++){
              idx3(w) = w;
            }
            arma::vec idx4(count2);
            for (int h = 0; h < count2; h++){
              idx4(h) = idx3(idx2(h));
            }
            arma::vec A2tmp(count2);
            for (int r = 0; r < count2; r++){
              A2tmp(r) = A2(idx4(r));
            }
            M(k, s, j) += sum(A2tmp * Phi(t, k) * Phi(t, s));
          }
          M(s, k, j) = M(k, s, j);
        }
      }
      M.slice(j) = M.slice(j) / sigma(j);
    }
  }

  Rcpp::List Mj_list;
  Mj_list["M"] = M;
  Mj_list["h"] = h;

  return(Mj_list);
}



// [[Rcpp::export]]

arma::mat V_solver(int K, int num_features, arma::mat self_V, arma::mat h,
                   arma::cube M) {
  for (int k = 0; k < K; k++) {
    arma::mat Mk(K, num_features);
    arma::vec a(num_features);
    arma::vec b(num_features);
    Mk = M.row(k);
    a = trans(Mk.row(k));
    b.zeros(num_features);
    b = trans(h.row(k));
    for (int s = 0; s < K; s++) {
      if (s != k) {
        b -= trans(Mk.row(s)) % self_V.col(s);
      }
    }
    double laplace_lam;
    laplace_lam = binary_search(num_features, 1, a, b);
    arma::vec lam(num_features);
    lam.fill(laplace_lam);
    self_V.col(k) = b / (a + lam);
  }
  return self_V;
}



// [[Rcpp::export]]

Rcpp::List muZ_transform(int dimZ1, int k, int selfK, arma::vec mutrans,
                         arma::mat Ztrans, arma::mat vec, arma::mat Z0,
                         arma::mat tmp1, arma::cube mat) {
  k -= 1;
  arma::mat mattmp(selfK, selfK);
  arma::vec tmp2(selfK);
  arma::vec tmp3(selfK - 1);
  for (int i = 0; i < dimZ1; i++) {
    mattmp = mat.row(i);
    tmp2 = trans(mattmp.row(k));
    int idx = 0;
    for (int j = 0; j < selfK; j++) {
      if (j != k){
        tmp3(idx) = tmp2(j);
        idx++;
      }
    }
    mutrans(i) = (vec(i, k) - sum(tmp3 % trans(tmp1.row(i)))) /
      sqrt(mat(i, k, k));
    Ztrans.row(i) = Z0.row(i) * sqrt(mat(i, k, k));
  }

  Rcpp::List trans_list;
  trans_list["mutrans"] = mutrans;
  trans_list["Ztrans"] = Ztrans;

  return(trans_list);
}



// [[Rcpp::export]]

arma::mat kronecker_muPhi(int K, int num_subjects, int num_times,
                          arma::mat Phi, arma::mat mu, arma::mat A0) {
  for (int i = 0; i < K; i++) {
    A0.col(i) = arma::kron(mu.col(i), Phi.col(i));
  }
  return(A0);
}



// [[Rcpp::export]]

arma::mat sigma_helper(int num_subjects, int num_times, int K, int idx1length,
                       arma::vec idx1, arma::mat A0, arma::mat A1,
                       arma::cube cov, arma::mat O1, arma::mat Phi) {
  arma::vec x(num_subjects);
  arma::mat xtmp1(num_subjects, K);
  for (int s = 0; s < K; s++) {
    for (int t = 0; t < K; t++) {
      if (t <= s) {
        arma::vec tmp1(idx1length);
        arma::vec tmp2(idx1length);
        for (int i = 0; i < idx1length; i++) {
          tmp1(i) = A0(idx1(i), s);
          tmp2(i) = A0(idx1(i), t);
        }
        A1(s, t) = sum(tmp1 % tmp2);
        xtmp1 = cov.slice(t);
        x = xtmp1.col(s);
        for (int j = 0; j < num_times; j++) {
          int count = 0;
          for (int m = 0; m < num_subjects; m++) {
            if (O1(m, j) == 1) {
              count++;
            }
          }
          arma::vec idx2(count);
          int index = 0;
          for (int n = 0; n < num_subjects; n++) {
            if (O1(n, j) == 1) {
              idx2(index) = n;
              index++;
            }
          }
          arma::vec idx3(num_subjects);
          arma::vec idx4(count);
          for (int l = 0; l < num_subjects; l++) {
            idx3(l) = l;
          }
          for (int k = 0; k < count; k++) {
            idx4(k) = idx3(idx2(k));
          }
          arma::vec xtmp2(count);
          for (int q = 0; q < count; q++) {
            xtmp2(q) = x(idx4(q));
          }
          A1(s, t) += sum(xtmp2 * Phi(j, s) * Phi(j, t));
        }
        A1(t, s) = A1(s, t);
      }
    }
  }
  return(A1);
}



// [[Rcpp::export]]

arma::vec sigma_solver(int num_subjects, int num_times, int num_features,
                       int idxlength, arma::mat XJ, arma::mat A0, arma::mat A1,
                       arma::mat V, arma::vec idx) {
  arma::vec x1(num_subjects * num_times);
  arma::vec x2(num_subjects * num_times);
  arma::vec x1tmp(idxlength);
  arma::vec x2tmp(idxlength);
  arma::vec sigma(num_features);
  double tmp;
  double num = 0.0;
  double den = 0.0;
  for (int i = 0; i < num_features; i++) {
    x1 = trans(XJ.row(i));
    x2 = A0 * trans(V.row(i));
    tmp = sum((V.row(i) * A1) % V.row(i));
    for (int j = 0; j < idxlength; j++) {
      x1tmp(j) = x1(idx(j));
      x2tmp(j) = x2(idx(j));
    }
    num += sum(square(x1tmp)) - 2 * sum(x1tmp % x2tmp) + tmp;
    den += idxlength;
  }
  for (int s = 0; s < num_features; s++) {
    sigma(s) = num/den;
  }
  return(sigma);
}



// [[Rcpp::export]]

arma::mat sigmaF_helper(int num_subjects, int K, bool fit_intercept,
                        arma::vec intercepts, arma::mat Z, arma::mat beta) {
  arma::mat fit(num_subjects, K);
  arma::mat fittmp(num_subjects, K);
  fittmp = Z * beta;
  if (fit_intercept) {
    for (int i = 0; i < K; i++) {
      fit.col(i) = fittmp.col(i) + intercepts(i);
    }
  }
  return(fit);
}



// [[Rcpp::export]]

arma::vec sigmaF_solver(int num_subjects, int K, arma::mat fit, arma::mat mu,
                        arma::cube cov) {
  arma::vec sigmaF(K);
  arma::vec tmp(K);
  arma::vec covdiag(K);
  arma::mat covtmp(K, K);
  for (int i = 0; i < num_subjects; i++) {
    tmp = trans(mu.row(i) - fit.row(i));
    covtmp = cov.row(i);
    for (int j = 0; j < K; j++) {
      covdiag(j) = covtmp(j, j);
    }
    sigmaF += square(tmp) + covdiag;
  }
  return(sigmaF);
}



// [[Rcpp::export]]

arma::mat left_creator(int h, arma::vec T1, arma::mat left) {
  double a;
  for (int t = 0; t < h - 1; t++) {
    a = 1.0 / (T1(t+1) - T1(t));
    left(t, t) = a;
    left(t, t+1) = -a;
  }
  return left;
}



// [[Rcpp::export]]

Rcpp::List ts_creator(int num_subjects, int num_times, arma::vec T1,
                      arma::mat O1) {
  arma::mat t(num_subjects, num_times);
  t.zeros(num_subjects, num_times);
  arma::mat s(num_subjects, num_times);
  s.zeros(num_subjects, num_times);
  for (int i = 0; i < num_subjects; i++) {
    arma::vec t0(num_times);
    for (int j = 0; j < num_times; j++) {
      t0(j) = j + 1;
    }
    t.row(i) = trans(t0);
  }
  for (int i = 0; i < num_times; i++) {
    arma::vec s0(num_subjects);
    for (int j = 0; j < num_subjects; j++) {
      s0(j) = j + 1;
    }
    s.col(i) = s0;
  }

  int count = 0;
  for (int i = 0; i < num_subjects; i++) {
    for (int j = 0; j < num_times; j++) {
      if (O1(i, j) == 1) {
        count++;
      }
    }
  }

  arma::vec ts(count);
  arma::vec s0(count);
  int idx = 0;
  for (int i = 0; i < num_subjects; i++) {
    for (int j = 0; j < num_times; j++) {
      if (O1(i, j) == 1) {
        ts(idx) = t(i, j);
        s0(idx) = s(i, j);
        idx++;
      }
    }
  }

  Rcpp::List ts_list;
  ts_list["ts"] = ts;
  ts_list["s0"] = s0;
  return ts_list;
}



// [[Rcpp::export]]

Rcpp::List mean_curve_update(int num_subjects, int num_times, int num_features,
                             double lam, arma::mat Psi, arma::mat O1,
                             arma::mat Omega, arma::cube X) {
  arma::mat B0(num_times, num_features);
  B0.zeros(num_times, num_features);
  arma::mat A0(num_times, num_times);
  A0 = trans(Psi) * Psi;
  int count = 0;
  for (int i = 0; i < num_subjects; i++) {
    for (int t = 0; t < num_times; t++) {
      if (O1(i, t) == 1) {
        count++;
      }
    }
  }

  arma::mat X0(count, num_features);
  X0.zeros(count, num_features);
  arma::mat A(num_times, num_times);
  arma::vec B(num_times);
  for (int j = 0; j < num_features; j++) {
    arma::mat x1(num_subjects, num_times);
    arma::vec x2(count);
    int idx = 0;
    x1 = X.slice(j);
    for (int i = 0; i < num_subjects; i++) {
      for (int t = 0; t < num_times; t++) {
        if (O1(i, t) == 1) {
          x2(idx) = x1(i, t);
          idx++;
        }
      }
    }
    X0.col(j) = x2;
    A = A0 + lam * Omega;
    B = trans(Psi) * x2;
    B0.col(j) = inv(A) * B;
  }

  Rcpp::List mean_curve;
  mean_curve["B0"] = B0;
  mean_curve["X0"] = X0;
  return mean_curve;
}



// [[Rcpp::export]]

arma::cube cv_creator(int nlam, int num_subjects, int num_times,
                      int num_features, arma::vec lams, arma::vec s0,
                      arma::mat O1, arma::mat Psi, arma::mat Omega,
                      arma::cube X, arma::cube cv) {

  for (int n = 0; n < nlam; n++) {
    double lam;
    lam = lams(n);
    arma::mat B0(num_times, num_features);
    B0.zeros(num_times, num_features);
    arma::mat A0(num_times, num_times);
    A0 = trans(Psi) * Psi;
    int count = 0;
    for (int i = 0; i < num_subjects; i++) {
      for (int t = 0; t < num_times; t++) {
        if (O1(i, t) == 1) {
          count++;
        }
      }
    }

    arma::mat X0(count, num_features);
    X0.zeros(count, num_features);
    arma::mat A(num_times, num_times);
    arma::vec B(num_times);
    for (int j = 0; j < num_features; j++) {
      arma::mat x1(num_subjects, num_times);
      arma::vec x2(count);
      int idx = 0;
      x1 = X.slice(j);
      for (int i = 0; i < num_subjects; i++) {
        for (int t = 0; t < num_times; t++) {
          if (O1(i, t) == 1) {
            x2(idx) = x1(i, t);
            idx++;
          }
        }
      }
      X0.col(j) = x2;
      A = A0 + lam * Omega;
      B = trans(Psi) * x2;
      B0.col(j) = inv(A) * B;
    }

    arma::mat R(count, num_features);
    arma::vec rj(count);
    R = X0 - Psi * B0;
    A0 = trans(Psi) * Psi + lam * Omega;
    A0 = inv(A0);
    for (int j = 0; j < num_features; j++) {
      rj = R.col(j);
      for (int i = 0; i < num_subjects; i++) {
        int s_count = 0;
        for (int m = 0; m < count; m++) {
          if (s0(m) == i) {
            s_count++;
          }
        }
        arma::vec si(s_count);
        int idx1 = 0;
        for (int z = 0; z < count; z++) {
          if (s0(z) == i) {
            si(idx1) = z;
            idx1++;
          }
        }
        arma::mat Psii(s_count, num_times);
        arma::vec tmp(num_times);
        for (int l = 0; l < s_count; l++) {
          int idx2;
          idx2 = si(l);
          tmp = trans(Psi.row(idx2));
          Psii.row(l) = trans(tmp);
        }
        arma::mat Vi(s_count, s_count);
        arma::mat diagV(s_count, s_count);
        arma::vec rij(s_count);
        Vi = diagV.eye(s_count, s_count) - (Psii * A0) * trans(Psii);
        Vi = inv(Vi);
        for (int h = 0; h < s_count; h++) {
          int idx3;
          idx3 = si(h);
          rij(h) = rj(idx3);
        }
        cv(i, j, n) = 1.0 / s_count * sum(square(Vi * rij));
      }
    }
  }

  return cv;
}



// [[Rcpp::export]]

arma::cube R_creator(int num_subjects, int num_times, int num_features,
                     arma::mat basis, arma::mat B0, arma::cube R) {
  arma::mat tmp(num_times, num_features);
  tmp = basis * B0;
  for (int j = 0; j < num_features; j++) {
    for (int i = 0; i < num_subjects; i++) {
      arma::mat Rj(num_subjects, num_times);
      Rj = R.slice(j);
      Rj.row(i) -= trans(tmp.col(j));
      R.slice(j) = Rj;
    }
  }
  return R;
}
