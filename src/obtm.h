/**
 * OBTM model of online BTM
 * Author: Xiaohui Yan(xhcloud@gmail.com)
 * 2013-6-6
 */
#ifndef _OBTM_H
#define _OBTM_H

#include <Rcpp.h>
#include <vector>

#include "biterm.h"
#include "pvec.h"
#include "pmat.h"

using namespace std;

class OBTM {
public:
  int K;				// topic number
  int W;				// vocabulary size
  int n_iter;			// maximum number of iteration of Gibbs Sampling

  double lam;					// decay coefficient
  Pvec<double> alpha;	// hyperparameters of p(z)
  Pmat<double> beta;			// hyperparameters of p(w|z)
  // record the sum over W in beta to save computation, size Kx1
  Pvec<double> beta_sum;

  // sample recorders
  Pvec<int> nb_z;   	// n(b|z), size K*1
  Pmat<int> nwz;	    // n(w,z), size M*K

  vector<Biterm> bs;  // training biterms
  
  bool has_background; 
  
  Pvec<double> pw_b;   // the background word distribution  

public:
  OBTM(int K, int W, double a, double b, int n_iter, double l, bool has_b = false):
  K(K), W(W), n_iter(n_iter), lam(l), has_background(has_b) {
    alpha.resize(K, a);
    beta.resize(K, W, b);
    pw_b.resize(W);
    beta_sum = beta.rowSum();
    nb_z.resize(K);
    nwz.resize(K, W);
  };

  void run(string input_dir, int n_day, string res_dir);

// private:
  void model_init();
  Rcpp::NumericVector split(int x, int n);
  Rcpp::NumericVector to_nvec(Pvec<double>  s);
  Rcpp::NumericMatrix to_nmat(Pmat<double>  s);
  void load_docs(string dfile);
  void proc_part(string pt);
  double loglik();

  void prepare_part();

  // run a batch Gibbs sampling precedure
  void Gibbs_sampling();

  // update estimate of a biterm
  void update_biterm(Biterm& bi);

  // reset topic proportions for biterm b
  void reset_biterm_topic(Biterm& bi);

  // assign topic proportions for biterm b
  void assign_biterm_topic(Biterm& bi, int k);

  // compute condition distribution p(z|b)
  void compute_pz_b(Biterm& bi, Pvec<double>& p);

  void save_res(string dir);
  void save_pz(string pt);
  void save_pw_z(string pt);
};

#endif
