#include <iostream>
#include <fstream>
#include <Rcpp.h>

#include "sampler.h"
#include "str_util.h"

#include "doc.h"
#include "ibtm.h"


// Function that prints
// the required sequence
Rcpp::NumericVector IBTM::split(int x, int n) 
{
  // int* arr = new int[n];
  Rcpp::NumericVector arr(n);
  // If we cannot split the
  // number into exactly 'N' parts
  if(x < n) {
    //cout<<"-1"<<" ";
    arr[1] = x;
  }
  // If x % n == 0 then the minimum
  // difference is 0 and all
  // numbers are x / n
  else if (x % n == 0)
  {
    for(int i=0;i<n;i++) {
      //cout<<(x/n)<<" ";
      arr[i] = x/n;
    }
  }
  else
  {
    
    // upto n-(x % n) the values
    // will be x / n
    // after that the values
    // will be x / n + 1
    int zp = n - (x % n);
    int pp = x/n;
    for(int i=0;i<n;i++)
    {
      
      if(i>= zp) {
        //cout<<(pp + 1)<<" ";
        arr[i] = pp+1;
      }
      else{
        //cout<<pp<<" ";
        arr[i] = pp;
      }
    }
  }
  return arr;
}


/**
 *   input_dir    contains docs orangnized by days start from 0,
 *                as day0.doc_wids, day1.day_wids, day2.day_wids, ...
 *   n_day        the number of days (starting from 0)
 *   output_dir   stores the p(z) and p(w|z) results for each day
 */
void IBTM::run(string input_dir, int n_day, string res_dir) {
  for (int d = 0; d < n_day; ++d) {
	Rcpp::Rcout << "----- proc day " << d << " ------" << endl;
	string pt = input_dir + str_util::itos(d) + ".txt";
	proc_part(pt);
	string dir = res_dir + "k" + str_util::itos(K) + ".day" + str_util::itos(d) + ".";
	save_res(dir);
  }
}


void IBTM::proc_part(string pt) {
  Rcpp::Rcout << "load docs: " << pt << endl;
  ifstream rf(pt.c_str());
  if (!rf) {
    Rcpp::Rcout << "file not find:" << pt << endl;
    Rcpp::stop(pt);;
  }

  string line;
  int n_d = 0;
  while(getline(rf, line)) {
	if (n_d % 1000 == 0) {
	  Rcpp::Rcout << "\rproc doc:" << n_d;
	}

	Doc doc(line);
	vector<Biterm> bis;
	doc.gen_biterms(bis);

    for (int b = 0; b < bis.size(); ++b)
      proc_biterm(bis[b]);

	++n_d;
  }
}

void IBTM::proc_biterm(Biterm& bi) {
  update_biterm(bi);

  // rejuvenate
  vector<int> idxs;
  gen_rej_idx(idxs);

  for (int i = 0; i < idxs.size(); ++i)
	update_biterm(bs[idxs[i]]);

  // add to rejuvenated list
  if (bs.size() < win) {
	bs.push_back(bi);
  }
  else {
	assert( bs.size() == win );
	bs[n_b % win] = bi;
  }

  // update the biterm counter
  ++n_b;
}

//select biterms index from [0, L-1] for rejevantation
// using fixing window + ... distribution
void IBTM::gen_rej_idx(vector<int>& idxs) {
  // for samll L, less rejuvenate
  int L = bs.size();
  for (int n = 0; n < min(L, n_rej); ++n) {
	int idx = Sampler::uni_sample(L);
	idxs.push_back(idx);
  }
}

// sample procedure for ith biterm
void IBTM::update_biterm(Biterm& bi) {
  if (bi.get_z() != -1)
	reset_biterm_topic(bi);

  // compute p(z|b)
  Pvec<double> pz;
  compute_pz_b(bi, pz);

  // sample topic for biterm b
  int k = Sampler::mult_sample(pz.to_vector());
  assign_biterm_topic(bi, k);
}

// reset topic assignment of biterm i
void IBTM::reset_biterm_topic(Biterm& bi) {
  int w1 = bi.get_wi();
  int w2 = bi.get_wj();
  int k = bi.get_z();

  --nb_z[k];				// update proportion of biterms in topic K
  --nwz[k][w1];			// update proportion w1's occurrence times in topic K
  --nwz[k][w2];
  assert(nb_z[k]>=0 && nwz[k][w1]>=0 && nwz[k][w2]>=0);

  bi.reset_z();
}

// compute p(z|w_i, w_j)
void IBTM::compute_pz_b(Biterm& bi, Pvec<double>& pz) {
  pz.resize(K);

  int w1 = bi.get_wi();
  int w2 = bi.get_wj();

  double pw1k, pw2k;
  for (int k = 0; k < K; ++k) {
    if (has_background && k == 0) {
      pw1k = pw_b[w1];
      pw2k = pw_b[w2];
    }
    else {
      pw1k = (nwz[k][w1] + beta) / (2 * nb_z[k] + W * beta);
      pw2k = (nwz[k][w2] + beta) / (2 * nb_z[k] + 1 + W * beta);
    }
	pz[k] = (nb_z[k] + alpha) * pw1k * pw2k;
  }
  pz.normalize();
}

// assign topic k to biterm i
void IBTM::assign_biterm_topic(Biterm& bi, int k) {
  int w1 = bi.get_wi();
  int w2 = bi.get_wj();
  bi.set_z(k);

  ++nb_z[k];
  ++nwz[k][w1];
  ++nwz[k][w2];
}

void IBTM::save_res(string dir) {
  string pt = dir + "pz";
  Rcpp::Rcout << "\nwrite p(z): " << pt << endl;
  save_pz(pt);

  string pt2 = dir + "pw_z";
  Rcpp::Rcout << "write p(w|z): " << pt2 << endl;
  save_pw_z(pt2);
}

// p(z) is determinated by the overall proportions
// of biterms in it
void IBTM::save_pz(string pt) {
  Pvec<double> pz(nb_z);
  pz.normalize(alpha);
  pz.write(pt);
}

void IBTM::save_pw_z(string pt) {
  Pmat<double> pw_z(nwz);
  pw_z.normr(beta);
  pw_z.write(pt);
}

double IBTM::loglik() { 
  Pvec<double> pz;
  pz.resize(K);
  for (int k = 0; k < K; k++){
    pz[k] = (nb_z[k] + alpha);
  }
  pz.normalize();
  
  Pmat<double> pw_z;
  pw_z = nwz.toDouble() + beta;
  pw_z.normr();
  // Pvec<double> lik_k(bs.size());
  double lik = 0;
  for (int b = 0; b < bs.size(); ++b) {
    double lik_k = 0;
    int w1 = bs[b].get_wi();
    int w2 = bs[b].get_wj();
    
    for (int k = 0; k < K; ++k) {
      // lik_k[b] += pw_z[k][w1] * pw_z[k][w2] * pz[k] ;
      lik_k += (pw_z[k][w1] * pw_z[k][w2] * pz[k] );
    }
    lik += log(lik_k);
  }
  // for (int k = 0; k < K; ++k) {
  // for (int b = 0; b < bs.size(); ++b) {
  //  lik += lik_k[b] ;
  //}
  double loglik;
  loglik = lik;
  return loglik;
}


