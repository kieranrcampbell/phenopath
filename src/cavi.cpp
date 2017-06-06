// [[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <cmath>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>

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
NumericMatrix update_greek_sum(int g, int p,
                               NumericMatrix greek_sum, 
                               double old_greek,
                               double new_greek,
                               NumericMatrix x) {
  /* This function strips out the alpha/beta value
   * from the previous iteration for the CAVI update
   */
  int N = greek_sum.ncol();
  // NumericMatrix new_greek_sum(clone(greek_sum));
  
  for(int i = 0; i < N; i++) {
    greek_sum(g,i) += -old_greek * x(i,p) + new_greek * x(i,p);
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
        for(int pp = 0; pp < P; pp++) {
          if(pp != p)
            sqe(g, i) += m_g(p,g) * m_g(pp,g) * x(i,p) * x(i,pp);
        }
      }
    }
  }
  return sqe;
}

// [[Rcpp::export]]
double calculate_fg(int g, NumericMatrix y, // NumericMatrix x, 
                    NumericVector m_z, NumericVector s_z, 
                    NumericVector m_lambda, NumericVector s_lambda,
                    NumericVector m_mu, NumericVector s_mu,
                    NumericMatrix alpha_sum, NumericMatrix beta_sum,
                    NumericMatrix alpha_square_sum, NumericMatrix beta_square_sum) {
  int N = y.nrow();
  
  double fg = 0.0;
  for(int i = 0; i < N; i++) {
    fg += m_mu[g] * m_mu[g] + s_mu[g]; // 1
    fg += 2 * m_mu[g] * alpha_sum(g,i); // 2
    fg += 2 * m_mu[g] * m_z[i] * m_lambda[g]; // 3
    fg += 2 * m_mu[g] * m_z[i] * beta_sum(g,i); // 4
    fg -= 2 * y(i,g) * m_mu[g]; // 5
    fg += alpha_square_sum(g,i); // 6
    fg += 2 * m_z[i] * m_lambda[g] * alpha_sum(g,i); // 7
    fg += 2 * m_z[i] * alpha_sum(g,i) * beta_sum(g,i); // 8
    fg -= 2 * y(i,g) * alpha_sum(g,i); // 9
    fg += (m_z[i] * m_z[i] + s_z[i]) * (m_lambda[g] * m_lambda[g] + s_lambda[g]); // 10
    fg += 2 * (m_z[i] * m_z[i] + s_z[i]) * m_lambda[g] * beta_sum(g,i); // 11
    fg -= 2 * m_lambda[g] * m_z[i] * y(i,g); // 12
    fg += (m_z[i] * m_z[i] + s_z[i]) * beta_square_sum(g,i);// 13
    fg -= 2 * y(i,g) * m_z[i] *  beta_sum(g,i); // 14
    fg += y(i,g) * y(i,g); // 15
    
  }
  return(fg);
}

// [[Rcpp::export]]
NumericMatrix cavi_update_z(NumericMatrix y, NumericMatrix x, 
                                    NumericVector m_lambda, NumericVector m_mu,
                                    NumericVector s_lambda, NumericMatrix m_alpha, 
                                    NumericMatrix m_beta, NumericMatrix s_beta,
                                    NumericVector a_tau, NumericVector b_tau,
                                    NumericVector q, double tau_q) {
  /***
   * This function returns an N-by-2 matrix where the entry in the 
   * i^th row is the update values of m_z_i and s_z_i^2 respectively
   */
  
  int N = y.nrow();
  int G = y.ncol();

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
      double pst_update_ig = (m_lambda[g] + beta_sum(g,i));
      pst_update_ig *= (y(i,g) - m_mu[g] - alpha_sum(g,i));
      pst_update_ig *= a_tau[g] / b_tau[g];
      pst_update(i,0) += pst_update_ig;      
    }
    
    // Calculate denominator
    for(int g = 0; g < G; g++) {
      double s_tmp = m_lambda[g] * m_lambda[g] + s_lambda[g];
      s_tmp +=  ((2 * m_lambda[g] * beta_sum(g,i)) + beta_sq_exp(g,i));
      
      pst_update(i, 1) += a_tau[g] / b_tau[g] * s_tmp;
    }
    
    pst_update(i, 0) /= pst_update(i, 1);
    pst_update(i, 1) = 1 / pst_update(i, 1); // invert to get variance
  }
  
  return pst_update;
}

// [[Rcpp::export]]
NumericMatrix cavi_update_mu(NumericMatrix y, NumericMatrix x, 
                             NumericVector m_z, NumericVector m_lambda,
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
      mu_update(g, 0) += y(i,g) - alpha_sum(g,i) - m_z[i] * (m_lambda[g] + beta_sum(g,i));
    }
    mu_update(g, 0) *= a_tau[g] / b_tau[g];
    
    mu_update(g, 0) /= mu_update(g, 1);
    mu_update(g, 1) = 1 / mu_update(g, 1);
  }

  
  return(mu_update);
}


// [[Rcpp::export]]
NumericMatrix cavi_update_lambda(NumericMatrix y, NumericMatrix x, 
                            NumericVector m_z, NumericVector s_z,
                            NumericMatrix m_alpha, NumericMatrix m_beta, 
                            NumericVector a_tau, NumericVector b_tau,
                            NumericVector m_mu, double tau_c) {
  
  int N = y.nrow();
  int G = y.ncol();
  
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
  NumericMatrix c_update(G, 2);
  fill(c_update.begin(), c_update.end(), 0.0);
  
  // calculate s_lambda 
  NumericVector m_s_square(N);
  double m_s_square_sum = 0.0;
  for(int i = 0; i < N; i++) {
    m_s_square[i] = m_z[i] * m_z[i] + s_z[i];
    m_s_square_sum += m_s_square[i];
  }

  // std::cout << m_s_square_sum << std::endl;
  
  for(int g = 0; g < G; g++)
    c_update(g, 1) = 1 / (a_tau[g] / b_tau[g] * m_s_square_sum + tau_c);

  
  for(int g = 0; g < G; g++) {
    for(int i = 0; i < N; i++) {
      c_update(g, 0) += (m_z[i] * y(i,g) -
      m_z[i] * m_mu[g] - m_z[i] * alpha_sum(g,i) -
      m_s_square[i] * beta_sum(g, i));
    }
    c_update(g, 0) *= (a_tau[g] / b_tau[g]);
    c_update(g, 0) *= c_update(g, 1);
  }
  
  return c_update;  
}

// [[Rcpp::export]]
NumericMatrix cavi_update_tau(NumericMatrix y, NumericMatrix x, 
                              NumericVector m_z, NumericVector s_z,
                              NumericVector m_lambda, NumericVector s_lambda,
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
    
    double fg = calculate_fg(g, y, m_z,  s_z, 
                             m_lambda, s_lambda, m_mu, s_mu,
                             alpha_sum, beta_sum,
                             alpha_square_sum, beta_square_sum);
    tau_update(g,1) = b + 0.5 * fg;
  }
  
  return tau_update;
}

// [[Rcpp::export]]
NumericVector cavi_update_alpha(NumericMatrix beta_sum, int p, int g, NumericMatrix y, NumericMatrix x, 
                                NumericVector m_z, NumericVector m_lambda,
                                NumericMatrix m_alpha, NumericMatrix m_beta,
                                NumericVector a_tau, NumericVector b_tau,
                                NumericVector m_mu, double tau_alpha) {
  /**
   * For alpha, beta and chi we update slightly differently - only for a given variable,
   * indexed by p and g
   */
  
  int N = y.nrow();
  int P = x.ncol();
  
  // NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  
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
    m_alpha_pg += x(i,p) * y(i,g);
    m_alpha_pg -= x(i,p) * m_mu[g];
    m_alpha_pg -= x(i,p) * m_z[i] * m_lambda[g];
    m_alpha_pg -= x(i,p) * m_z[i] * beta_sum(g,i);
    m_alpha_pg -= x(i,p) * alpha_sum_no_p[i];
  }
  m_alpha_pg *= a_tau[g] / b_tau[g];
  
  return NumericVector::create(m_alpha_pg * s_alpha_pg, s_alpha_pg);
}


// [[Rcpp::export]]
NumericVector cavi_update_beta(NumericMatrix alpha_sum, int p, int g, NumericMatrix y, NumericMatrix x, 
                                NumericVector m_z, NumericVector s_z, NumericVector m_lambda,
                                NumericMatrix m_alpha, NumericMatrix m_beta,
                                NumericVector a_tau, NumericVector b_tau,
                                NumericMatrix a_chi, NumericMatrix b_chi,
                                NumericVector m_mu) {

  int N = y.nrow();
  int P = x.ncol();
  
  // NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  
  /** 
   * We start by calculating some useful quantities
   */
  NumericVector ms_vec(N);
  for(int i = 0; i < N; i++) {
    ms_vec[i] = m_z[i] * m_z[i] + s_z[i];
  }
  
  NumericVector beta_sum_no_p(N, 0.0);
  for(int i = 0; i < N; i++) {
    for(int pp = 0; pp < P; pp++) {
      if(pp != p)
        beta_sum_no_p[i] += m_beta(pp,g) * x(i,pp);
    }
  }


  // Calculate s_beta_pg
  double s_beta_pg = 0.0;
  for(int i = 0; i < N; i++) {
    s_beta_pg += (ms_vec[i] * x(i,p) * x(i,p));
  }
  s_beta_pg *= a_tau[g] / b_tau[g];
  s_beta_pg += a_chi(p,g) / b_chi(p,g);
  s_beta_pg = 1 / s_beta_pg;
  
  double m_beta_pg = 0.0;
  
  for(int i = 0; i < N; i++) {
    m_beta_pg += x(i,p) * m_z[i] * y(i,g);
    m_beta_pg -= x(i,p) * m_lambda[g] * ms_vec[i];
    m_beta_pg -= m_z[i] * alpha_sum(g,i) * x(i,p);
    m_beta_pg -= ms_vec[i] * beta_sum_no_p[i] * x(i,p);
  }
  m_beta_pg *= (a_tau[g] / b_tau[g]) * s_beta_pg;
  return NumericVector::create(m_beta_pg, s_beta_pg);
  
}


// [[Rcpp::export]]
NumericVector cavi_update_chi(double m_beta_pg, double s_beta_pg,
                              double a_beta, double b_beta) {
  double a_new = a_beta + 0.5;
  double b_new = b_beta + 0.5 * (m_beta_pg * m_beta_pg + s_beta_pg);
  
  return NumericVector::create(a_new, b_new);
}

// [[Rcpp::export]]
double calculate_E_log_Y_given_theta(NumericMatrix y, NumericMatrix x, 
                                     NumericVector m_z, NumericVector s_z, 
                                     NumericVector m_lambda, NumericVector s_lambda,
                                     NumericMatrix m_alpha, NumericMatrix s_alpha,
                                     NumericMatrix m_beta, NumericMatrix s_beta,
                                     NumericVector a_tau, NumericVector b_tau,
                                     NumericVector m_mu, NumericVector s_mu) {
  int N = y.nrow();
  //int P = x.ncol();
  int G = y.ncol();
  
  double ely = 0.0; // Expectation of log p(Y|\theta)
  double pi = 3.14159265359;
  
  // temporary holders
  NumericMatrix alpha_sum = calculate_greek_sum(m_alpha, x);
  NumericMatrix beta_sum = calculate_greek_sum(m_beta, x);
  NumericMatrix alpha_square_sum = greek_square_exp(m_alpha, s_alpha, x);
  NumericMatrix beta_square_sum = greek_square_exp(m_beta, s_beta, x);
  
  for(int g = 0; g < G; g++) {
    double e_log_tau = (boost::math::digamma(a_tau[g]) - log(b_tau[g]));
    ely += N / 2 * (e_log_tau - log(2 * pi));
    double fg = calculate_fg(g, y,m_z,  s_z, 
                 m_lambda, s_lambda, m_mu, s_mu,
                 alpha_sum, beta_sum,
                 alpha_square_sum, beta_square_sum);
    ely -= a_tau[g] / (2 * b_tau[g]) * fg;
  }
  return ely;
}

// [[Rcpp::export]]
double calculate_E_log_p(NumericVector m_z, NumericVector s_z, 
                             NumericVector m_lambda, NumericVector s_lambda,
                             NumericMatrix m_alpha, NumericMatrix s_alpha,
                             NumericMatrix m_beta, NumericMatrix s_beta,
                             NumericVector a_tau, NumericVector b_tau,
                             NumericVector m_mu, NumericVector s_mu,
                             NumericMatrix a_chi, NumericMatrix b_chi,
                             NumericVector q, double tau_q, double tau_mu, double tau_c,
                             double a, double b, double tau_alpha,
                             double a_beta, double b_beta) {
  int N = m_z.size();
  int G = m_lambda.size();
  int P = a_chi.nrow();
  
  double elp = 0.0;
  double pi = 3.14159265359;
  
  for(int i = 0; i < N; i++) {
    elp += 0.5 * log(tau_q / (2*pi)) - (tau_q / 2) * (m_z[i] * m_z[i] + s_z[i] - 2 * m_z[i] * q[i] + q[i] * q[i]);
  }
  
  for(int g = 0; g < G; g++) {
    elp += 0.5 * log(tau_mu / (2*pi)) - tau_mu / 2 * (m_mu[g] * m_mu[g] + s_mu[g]); 
    elp += 0.5 * log(tau_c / (2*pi)) - tau_c / 2 * (m_lambda[g] * m_lambda[g] + s_lambda[g]);
    
    elp += (a - 1) * (boost::math::digamma(a_tau[g]) - log(b_tau[g])) -
      a_tau[g] / b_tau[g] * b + a * log(b) - boost::math::lgamma(a);
    
    for(int p = 0; p < P; p++) {
      elp += 0.5 * log(tau_alpha / (2*pi)) - tau_alpha / 2 * (m_alpha(p,g) * m_alpha(p,g) + s_alpha(p,g)); 
      
      elp += 0.5 * (boost::math::digamma(a_chi(p,g)) - log(b_chi(p,g))) - 
        0.5 * log(2 * pi) -
        a_chi(p,g) / (2 * b_chi(p,g)) * (m_beta(p,g) * m_beta(p,g) + s_beta(p,g));
      
      elp += (a_beta - 1) * (boost::math::digamma(a_chi(p,g)) - log(b_chi(p,g))) -
        a_chi(p,g) / b_chi(p,g) * b_beta + a_beta * log(b_beta) - boost::math::lgamma(a_beta);
    }
  }
  
  return elp;
}

// [[Rcpp::export]]

double calculate_E_log_q(NumericVector s_z, NumericVector s_lambda,
                  NumericMatrix s_alpha, NumericMatrix s_beta,
                  NumericVector a_tau, NumericVector b_tau,
                  NumericVector s_mu,
                  NumericMatrix a_chi, NumericMatrix b_chi,
                  int model_mu) {
  double elq = 0.0;
  double pi = 3.14159265359;
  
  int N = s_z.size();
  int G = s_lambda.size();
  int P = a_chi.nrow();
  
  for(int i = 0; i < N; i++)
    elq -= 0.5 * log(s_z[i]) - 0.5 *  log(2*pi);
  
  for(int g = 0; g < G; g++) {
    if(model_mu == 1) {
      elq -= 0.5 * log(s_mu[g]) - 0.5 *  log(2 * pi); 
    }
    elq -= 0.5 * log(s_lambda[g]) - 0.5 * log(2 * pi);
    
    elq += (a_tau[g] - 1) *
      (boost::math::digamma(a_tau[g]) - log(b_tau[g])) - a_tau[g] +
      a_tau[g] * log(b_tau[g]) - boost::math::lgamma(a_tau[g]);
    
    for(int p = 0; p < P; p++) {
      elq -= 0.5 * log(s_alpha(p,g)) - 0.5 * log(2 * pi); 
      
      elq -= 0.5 * log(s_beta(p,g)) - 0.5 * log(2 * pi);
      
      elq += (a_chi(p,g) - 1) *
        (boost::math::digamma(a_chi(p,g)) - log(b_chi(p,g))) - a_chi(p,g) +
        a_chi(p,g) * log(b_chi(p,g)) - boost::math::lgamma(a_chi(p,g));
    }
  }
  
  return elq;
}

// [[Rcpp::export]]
NumericVector calculate_elbo(NumericMatrix y, NumericMatrix x, 
                      NumericVector m_z, NumericVector s_z, 
                      NumericVector m_lambda, NumericVector s_lambda,
                      NumericMatrix m_alpha, NumericMatrix s_alpha,
                      NumericMatrix m_beta, NumericMatrix s_beta,
                      NumericVector a_tau, NumericVector b_tau,
                      NumericMatrix a_chi, NumericMatrix b_chi,
                      NumericVector m_mu, NumericVector s_mu,
                      NumericVector q, double tau_q, double tau_mu, double tau_c,
                      double a, double b, double tau_alpha,
                      double a_beta, double b_beta, int model_mu) {
  
  double ely = calculate_E_log_Y_given_theta( y,  x,
                                              m_z,  s_z,
                                              m_lambda,  s_lambda,
                                              m_alpha,  s_alpha,
                                              m_beta,  s_beta,
                                              a_tau,  b_tau,
                                              m_mu,  s_mu);
  
  double elp = calculate_E_log_p( m_z,  s_z,
                                  m_lambda,  s_lambda,
                                  m_alpha,  s_alpha,
                                  m_beta,  s_beta,
                                  a_tau,  b_tau,
                                  m_mu,  s_mu,
                                  a_chi,  b_chi,
                                  q,  tau_q,  tau_mu,  tau_c,
                                  a,  b,  tau_alpha,
                                  a_beta,  b_beta);

  double elq = calculate_E_log_q( s_z,  s_lambda,
                 s_alpha,  s_beta,
                 a_tau,  b_tau,
                 s_mu,
                 a_chi,  b_chi, model_mu);

  
  return NumericVector::create(ely, elp, elq);
}


  
  