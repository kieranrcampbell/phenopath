// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix calculate_greek_sum(NumericMatrix greek, NumericMatrix x) {
  /**
   * Despite the slightly odd name, this function takes either alpha or beta as 'greek'
   * and calculates the G-by-N matrix \sum_p greek_pg x_ip
   */
  int N = x.nrow();
  int P = x.ncol();
  int G = greek.ncol();
  
  NumericMatrix greek_sum(G, N);
  fill(greek_sum.begin(), greek_sum.end(), 0.0);
  
  for(int i = 0; i < N; i++) {
    for(int g = 0; g < G; g++) {
      for(int p = 0; p < P; p++) {
        greek_sum(g, i) += greek(p,g) * x(i,p);
      }
    }
  }
  return greek_sum;
}

// [[Rcpp::export]]
NumericMatrix greek_square_exp(NumericMatrix m_g, NumericMatrix s_g, NumericMatrix x) {
  /**
   * this calculates E[(\sum_p greek_{pg} x_{ip})^2] = 
   * \sum_p (m_g_pg^2 + s_g_pg^2) x_ip^2 + \sum_{p,p':p\neqp'} m_g_pg m_g_p'g x_ip x_ip'
   */
  int P = m_g.nrow();
  int G = m_g.ncol();
  int N = x.nrow();
  
  NumericMatrix sqe(G, N);
  fill(sqe.begin(), sqe.end(), 0.0);
  
  for(int g = 0; g < G; g++) {
    for(int i = 0; i < N; i++) {
      for(int p = 0; p < P; p++) {
        sqe(g, i) += (m_g(p,g) * m_g(p,g) + s_g(p,g)) * x(i,p) * x(i,p);
        for(int pp = 0; pp < P; pp++) 
          if(pp != p)
            sqe(g, i) += m_g(p,g) * m_g(pp,g) * x(i,p) * x(i,pp);
      }
    }
  }
  return sqe;
}

// [[Rcpp::export]]
NumericMatrix cavi_update_pst(NumericMatrix y, NumericMatrix x, 
                                    NumericVector m_c, NumericVector m_mu,
                                    NumericVector s_c, NumericMatrix m_alpha, 
                                    NumericMatrix m_beta, NumericMatrix s_beta,
                                    NumericVector a_tau, NumericVector b_tau,
                                    NumericVector q, double tau_q) {
  /***
   * This function returns an N-by-2 matrix where the entry in the 
   * i^th row is the update values of m_t_i and s_t_i^2 respectively
   */
  
  int N = y.nrow();
  int G = y.ncol();
  //int P = x.ncol();

  NumericMatrix pst_update(N, 2);
  
  fill(pst_update.begin(), pst_update.end(), 0.0);
  
  /**
   * Note that both m and s^2 require \sum_p m_beta_pg x_ip and m
   * further requires \sum_p m_alpha_pg, so we compute those separately first
   */
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  NumericMatrix beta_sq_exp = greek_square_exp(m_beta, s_beta, x);
  
  for(int i = 0; i < N; i++) {
    pst_update(i, 0) = tau_q * q[i];
    pst_update(i, 1) = tau_q;
    
    // Calculate numerator
    for(int g = 0; g < G; g++) {
      pst_update(i, 0) += a_tau[g] / b_tau[g] * (
        m_c[g] + beta_sum(g,i)
      ) * (y(i,g) - m_mu[g] - alpha_sum(g,i));
    }
    
    // Calculate denominator
    for(int g = 0; g < G; g++) {
      double s_tmp = pow(m_c[g], 2.0) + s_c[g];
      s_tmp += 2 * m_c[g] * beta_sum(g,i) + beta_sq_exp(g,i);
      
      pst_update(i, 1) += a_tau[g] / b_tau[g] * s_tmp;
    }
    
    pst_update(i, 0) /= pst_update(i, 1);
    pst_update(i, 1) = 1 / pst_update(i, 1); // invert to get variance
  }
  
  return pst_update;
}

// [[Rcpp::export]]
NumericMatrix cavi_update_mu(NumericMatrix y, NumericMatrix x, 
                             NumericVector m_t, NumericVector m_c,
                             NumericMatrix m_alpha, NumericMatrix m_beta, 
                             NumericVector a_tau, NumericVector b_tau,
                             double tau_mu) {
  
  int N = y.nrow();
  int G = y.ncol(); 
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  NumericMatrix mu_update(G, 2);
  fill(mu_update.begin(), mu_update.end(), 0.0);
  
  // update s_mu
  for(int g = 0; g < G; g++)
    mu_update(g, 1) = a_tau[g] / b_tau[g] * N + tau_mu;
  
  // update m_mu
  for(int g = 0; g < G; g++) {
    for(int i = 0; i < N; i++) {
      mu_update(g, 0) += y(i,g) - alpha_sum(g,i) - m_t[i] * (m_c[g] + beta_sum(g,i));
    }
    mu_update(g, 0) *= a_tau[g] / b_tau[g];
    
    mu_update(g, 0) /= mu_update(g, 1);
    mu_update(g, 1) = 1 / mu_update(g, 1);
  }

  
  return(mu_update);
}


// [[Rcpp::export]]
NumericMatrix cavi_update_c(NumericMatrix y, NumericMatrix x, 
                            NumericVector m_t, NumericVector s_t,
                            NumericMatrix m_alpha, NumericMatrix m_beta, 
                            NumericVector a_tau, NumericVector b_tau,
                            NumericVector m_mu, double tau_c) {
  
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  // 
  NumericMatrix c_update(G, 2);
  fill(c_update.begin(), c_update.end(), 0.0);
  
  // calculate s_c 
  NumericVector m_s_square(N);
  double m_s_square_sum = 0.0;
  for(int i = 0; i < N; i++) {
    m_s_square[i] = pow(m_t[i], 2) + s_t[i];
    m_s_square_sum += m_s_square[i];
  }

  for(int g = 0; g < G; g++)
    c_update(g, 1) = 1 / (a_tau[g] / b_tau[g] * m_s_square_sum + tau_c);

  
  for(int g = 0; g < G; g++) {
    for(int i = 0; i < N; i++) {
      c_update(g, 0) += (m_t[i] * y(i,g)
      - m_t[i] * m_mu[g] - m_t[i] * alpha_sum(g,i) -
      m_s_square[i] * beta_sum(g, i));
    }
    c_update(g, 0) *= a_tau[g] / b_tau[g];
    c_update(g, 0) *= c_update(g, 1);
  }
  
  return c_update;  
}

// [[Rcpp::export]]
NumericMatrix cavi_update_tau(NumericMatrix y, NumericMatrix x, 
                              NumericVector m_t, NumericVector s_t,
                              NumericVector m_c, NumericVector s_c,
                              NumericMatrix m_alpha, NumericMatrix m_beta,
                              NumericMatrix s_alpha, NumericMatrix s_beta,
                              NumericVector m_mu, NumericVector s_mu,
                              double a, double b) {
  
  /***
   * Here we return a G-by-2 matrix, where the first column is the value of
   * a_tau and the second b_tau
   */
  
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  NumericMatrix alpha_square_sum = greek_square_exp(m_alpha, s_alpha, x);
  NumericMatrix beta_square_sum = greek_square_exp(m_beta, s_beta, x);
  
  NumericMatrix tau_update(G, 2);
  
  for(int g = 0; g < G; g++)
    tau_update(g,0) = a + N / 2.0;
  
  for(int g = 0; g < G; g++) {
    double fg = 0.0;
    for(int i = 0; i < N; i++) {
      fg += y(i,g) * y(i,g) + m_mu[g] * m_mu[g] + s_mu[g] + 
        (m_t[i] * m_t[i] + s_t[i]) * (m_c[g] * m_c[g] + s_c[g]);
      fg += alpha_square_sum(g,i) + (m_t[i] * m_t[i] + s_t[i]) * beta_square_sum(g,i);
      fg -= 2 * y(i,g) * (m_mu[g] + alpha_sum(g,i) + m_t[i] * (m_c[g] + beta_sum(g,i)));
      fg += 2 * (m_mu[g] * alpha_sum(g,i) + m_mu[g] * m_t[i] * m_c[g]);
      fg += 2 * (m_mu[g] * m_t[i] * beta_sum(g,i) + m_t[i] * m_c[g] * alpha_sum(g,i));
      fg += 2 * (m_t[i] * alpha_sum(g,i) * beta_sum(g,i) + 
        (m_t[i] * m_t[i] + s_t[i]) * m_c[g]*  beta_sum(g,i)  );
      }

    tau_update(g,1) = b + 0.5 * fg;
  }
  
  return tau_update;
}

// [[Rcpp::export]]
NumericVector cavi_update_alpha(int p, int g, NumericMatrix y, NumericMatrix x, 
                                NumericVector m_t, NumericVector m_c,
                                NumericMatrix m_alpha, NumericMatrix m_beta,
                                NumericVector a_tau, NumericVector b_tau,
                                NumericVector m_mu, double tau_alpha) {
  /**
   * For alpha, beta and chi we update slightly differently - only for a given variable,
   * indexed by p and g
   */
  
  int N = y.nrow();
  int P = x.ncol();
  
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  double s_alpha_pg = tau_alpha;
  for(int i = 0; i < N; i++)
    s_alpha_pg += a_tau[g] / b_tau[g] * x(i,p) * x(i,p);
  s_alpha_pg = 1.0 / s_alpha_pg;

  // need to calculate alpha sum without the p'th entry
  NumericVector alpha_sum_no_p(N, 0.0);
  for(int i = 0; i < N; i++) {
    for(int pp = 0; pp < P; pp++) {
      if(pp != p)
        alpha_sum_no_p[i] += m_alpha(pp,g) * x(i,pp);
    }
  }
  
  double m_alpha_pg = 0.0;
  for(int i = 0; i < N; i++) {
    m_alpha_pg += x(i,p) * (y(i,g) - m_mu[g] - m_t[i] * (m_c[g] + beta_sum(g,i)) - alpha_sum_no_p[i]);
  }
  m_alpha_pg *= a_tau[g] / b_tau[g];
  
  return NumericVector::create(m_alpha_pg * s_alpha_pg, s_alpha_pg);
}


// [[Rcpp::export]]
NumericVector cavi_update_beta(int p, int g, NumericMatrix y, NumericMatrix x, 
                                NumericVector m_t, NumericVector s_t, NumericVector m_c,
                                NumericMatrix m_alpha, NumericMatrix m_beta,
                                NumericVector a_tau, NumericVector b_tau,
                                NumericMatrix a_chi, NumericMatrix b_chi,
                                NumericVector m_mu) {

  int N = y.nrow();
  int P = x.ncol();
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  
  /** 
   * We start by calculating some useful quantities
   */
  NumericVector ms_vec(N);
  for(int i = 0; i < N; i++)
    ms_vec[i] = pow(m_t[i], 2) + s_t[i];
  
  NumericVector beta_sum_no_p(N, 0.0);
  for(int i = 0; i < N; i++) {
    for(int pp = 0; pp < P; pp++) {
      if(pp != p)
        beta_sum_no_p[i] += m_beta(pp,g) * x(i,pp);
    }
  }
  
  // Calculate s_beta_pg
  double s_beta_pg = a_chi(p,g) / b_chi(p,g);
  for(int i = 0; i < N; i++) {
    s_beta_pg += a_tau[g] / b_tau[g] * ms_vec[i] * pow(x(i,p), 2);
  }
  s_beta_pg = 1 / s_beta_pg;
  
  double m_beta_pg = 0.0;
  
  for(int i = 0; i < N; i++) {
    m_beta_pg += x(i,p) * ( 
      m_t[i] * y(i,g) - m_t[i] * m_mu[g] - ms_vec[i] * m_c[g] - m_t[i] * alpha_sum(g,i) - 
        ms_vec[i] * beta_sum_no_p[i]
    );
  }
  m_beta_pg *= a_tau[g] / b_tau[g] * s_beta_pg;
  return NumericVector::create(m_beta_pg, s_beta_pg);
  
}


// [[Rcpp::export]]
NumericVector cavi_update_chi(double m_beta_pg, double s_beta_pg,
                              double a_beta, double b_beta) {
  double a_new = a_beta + 0.5;
  double b_new = b_beta + 0.5 * (pow(m_beta_pg, 2) + s_beta_pg);
  
  return NumericVector::create(a_new, b_new);
}

// [[Rcpp::export]]
double calculate_E_log_Y_given_theta(NumericMatrix y, NumericMatrix x, 
                                     NumericVector m_t, NumericVector s_t, 
                                     NumericVector m_c, NumericVector s_c,
                                     NumericMatrix m_alpha, NumericMatrix s_alpha,
                                     NumericMatrix m_beta, NumericMatrix s_beta,
                                     NumericVector a_tau, NumericVector b_tau,
                                     NumericVector m_mu, NumericVector s_mu) {
  int N = y.nrow();
  //int P = x.ncol();
  int G = y.ncol();
  
  double ely = 0.0; // Expectation of log p(Y|\theta)
  
  // temporary holders
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  NumericMatrix alpha_square_sum = greek_square_exp(m_alpha, s_alpha, x);
  NumericMatrix beta_square_sum = greek_square_exp(m_beta, s_beta, x);
  
  for(int g = 0; g < G; g++) {
    ely += N / 2 * (boost::math::digamma(a_tau[g]) - log(b_tau[g]));
    double fg = 0.0;
    for(int i = 0; i < N; i++) {
      fg += y(i,g) * y(i,g) + m_mu[g] * m_mu[g] + s_mu[g] + 
        (m_t[i] * m_t[i] + s_t[i]) * (m_c[g] * m_c[g] + s_c[g]);
      fg += alpha_square_sum(g,i) + (m_t[i] * m_t[i] + s_t[i]) * beta_square_sum(g,i);
      fg -= 2 * y(i,g) * (m_mu[g] + alpha_sum(g,i) + m_t[i] * (m_c[g] + beta_sum(g,i)));
      fg += 2 * (m_mu[g] * alpha_sum(g,i) + m_mu[g] * m_t[i] * m_c[g]);
      fg += 2 * (m_mu[g] * m_t[i] * beta_sum(g,i) + m_t[i] * m_c[g] * alpha_sum(g,i));
      fg += 2 * (m_t[i] * alpha_sum(g,i) * beta_sum(g,i) + 
        (m_t[i] * m_t[i] + s_t[i]) * m_c[g]*  beta_sum(g,i)  );
    }
    ely -= a_tau[g] / 2 * b_tau[g] * fg;
  }
  return ely;
}

// [[Rcpp::export]]
double calculate_E_log_p(NumericVector m_t, NumericVector s_t, 
                             NumericVector m_c, NumericVector s_c,
                             NumericMatrix m_alpha, NumericMatrix s_alpha,
                             NumericMatrix m_beta, NumericMatrix s_beta,
                             NumericVector a_tau, NumericVector b_tau,
                             NumericVector m_mu, NumericVector s_mu,
                             NumericMatrix a_chi, NumericMatrix b_chi,
                             NumericVector q, double tau_q, double tau_mu, double tau_c,
                             double a, double b, double tau_alpha,
                             double a_beta, double b_beta) {
  int N = m_t.size();
  int G = m_c.size();
  int P = a_chi.nrow();
  
  double elp = 0.0;
  
  for(int i = 0; i < N; i++)
    elp -= tau_q / 2 * (m_t[i] * m_t[i] + s_t[i] - 2 * m_t[i] * q[i]);
  
  for(int g = 0; g < G; g++) {
    elp -= tau_mu / 2 * (m_mu[g] * m_mu[g] + s_mu[g]) + 
      tau_c / 2 * (m_c[g] * m_c[g] + s_c[g]);
    elp += (a - 1) * (boost::math::digamma(a_tau[g]) - log(b_tau[g])) -
      a_tau[g] / b_tau[g] * b;
    for(int p = 0; p < P; p++) {
      elp -= tau_alpha / 2 * (m_alpha(p,g) * m_alpha(p,g) + s_alpha(p,g)) +
        a_chi(p,g) / (2 * b_chi(p,g)) * (m_beta(p,g) * m_beta(p,g) + s_beta(p,g));
      elp += (a_beta - 1) * (boost::math::digamma(a_chi(p,g)) - log(b_chi(p,g))) -
        a_chi(p,g) / b_chi(p,g) * b_beta;
    }
  }
  
  return elp;
}

// [[Rcpp::export]]

double calculate_E_log_q(NumericVector s_t, NumericVector s_c,
                  NumericMatrix s_alpha, NumericMatrix s_beta,
                  NumericVector a_tau, NumericVector b_tau,
                  NumericVector s_mu,
                  NumericMatrix a_chi, NumericMatrix b_chi) {
  double elq = 0.0;
  
  int N = s_t.size();
  int G = s_c.size();
  int P = a_chi.nrow();
  
  for(int i = 0; i < N; i++)
    elq -= 0.5 * s_t[i];
  
  for(int g = 0; g < G; g++) {
    elq -= 0.5 * s_mu[g] + 0.5 * s_c[g] - (a_tau[g] - 1) * 
      (boost::math::digamma(a_tau[g]) - log(b_tau[g])) + a_tau[g];
    for(int p = 0; p < P; p++) {
      elq -= 0.5 * s_alpha(p,g) + 0.5 * s_beta(p,g) - (a_chi(p,g) - 1) * 
        (boost::math::digamma(a_chi(p,g)) - log(b_chi(p,g))) + a_chi(p,g);
    }
  }
  
  return elq;
}

// [[Rcpp::export]]
double calculate_elbo(NumericMatrix y, NumericMatrix x, 
                      NumericVector m_t, NumericVector s_t, 
                      NumericVector m_c, NumericVector s_c,
                      NumericMatrix m_alpha, NumericMatrix s_alpha,
                      NumericMatrix m_beta, NumericMatrix s_beta,
                      NumericVector a_tau, NumericVector b_tau,
                      NumericMatrix a_chi, NumericMatrix b_chi,
                      NumericVector m_mu, NumericVector s_mu,
                      NumericVector q, double tau_q, double tau_mu, double tau_c,
                      double a, double b, double tau_alpha,
                      double a_beta, double b_beta) {
  
  double ely = calculate_E_log_Y_given_theta( y,  x,
                                              m_t,  s_t,
                                              m_c,  s_c,
                                              m_alpha,  s_alpha,
                                              m_beta,  s_beta,
                                              a_tau,  b_tau,
                                              m_mu,  s_mu);
  
  double elp = calculate_E_log_p( m_t,  s_t,
                                  m_c,  s_c,
                                  m_alpha,  s_alpha,
                                  m_beta,  s_beta,
                                  a_tau,  b_tau,
                                  m_mu,  s_mu,
                                  a_chi,  b_chi,
                                  q,  tau_q,  tau_mu,  tau_c,
                                  a,  b,  tau_alpha,
                                  a_beta,  b_beta);

  double elq = calculate_E_log_q( s_t,  s_c,
                 s_alpha,  s_beta,
                 a_tau,  b_tau,
                 s_mu,
                 a_chi,  b_chi);
  
  return ely + elp - elq;
}












  
  