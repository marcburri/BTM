/**
 * IBTM algorithm of online BTM
 * Author: Xiaohui Yan(xhcloud@gmail.com)
 * 2013-6-6
 */
#ifndef _IBTM_H
#define _IBTM_H

#include <Rcpp.h>
#include <string>
#include <vector>

#include "biterm.h"
#include "pvec.h"
#include "pmat.h"

using namespace std;

class IBTM {
public:
  int K;				// topic number
  int W;				// vocabulary size
  int n_iter;			// maximum number of iteration of Gibbs Sampling

  double alpha;			// hyperparameters of p(z)
  double beta;			// hyperparameters of p(w|z)

  int win_nrej;		// size of the biterm sliding window
  int n_rej;	// number of biterms to rejuvenate
  long n_b;		// number of biterms processed

  // sample recorders
  Pvec<int> nb_z;	// n(b|z), size K*1
  Pmat<int> nwz;	  // n(w,z), size K*W

  vector<Biterm> bs;  // training biterms
  
  bool has_background; 
  
  Pvec<double> pw_b;   // the background word distribution  
  
public:
  IBTM(int K, int W, double a, double b, int win_nrej, int n_rej, bool has_b = false):
  K(K), W(W),  alpha(a), beta(b),
  win_nrej(win_nrej), n_rej(n_rej), n_b(0) , has_background(has_b) {  
    nb_z.resize(K);
    nwz.resize(K, W);
    pw_b.resize(W);
  }
  
  Rcpp::NumericVector split(int x, int n);
  
  void run(string input_dir, int n_day, string res_dir);
  
  void proc_part(string pt);
  
  void model_init();
  
// private:
  void proc_biterm(Biterm& bi);
  void gen_rej_idx(vector<int>& idxs);
  double loglik();

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
