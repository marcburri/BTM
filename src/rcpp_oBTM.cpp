#include <Rcpp.h>
#include <cstdlib>
#include <string.h>
#include <string>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include "sampler.h"


#include "obtm.h"
#include "infer.h"


using namespace std;

// [[Rcpp::export]]
SEXP obtm(Rcpp::List biterms, Rcpp::CharacterVector x, int K, int W, double a, double b, int iter, int win = 15, double lam = 1, int n_part = 10,int trace = 0, int check_convergence = 0, double convergence_tol = 0.001, bool background = false) {
  Rcpp::Function format_posixct("format.POSIXct");
  Rcpp::Function sys_time("Sys.time");
  int biterms_size = biterms.size();
  bool has_background = background;
  Rcpp::CharacterVector doc_ids = x.attr("names");
  std::string context_id;
  Rcpp::List context_id_biterm;
  Rcpp::IntegerVector term1;
  Rcpp::IntegerVector term2;
  Rcpp::IntegerVector cooc;

  Rcpp::XPtr<OBTM> obtm(new OBTM(K, W, a, b, iter, lam, has_background), true);
  std::string line;


  Rcpp::NumericVector arru(n_part);
  arru = obtm->split(x.size(),n_part);

  Rcpp::NumericVector arrl(n_part);

  int h = 0;
  for (int d = 0; d < n_part+1; ++d) {
    arrl[d] = h;
    h = arru[d] + 1;
    arru[d+1] += arru[d] ;
  }

  for (int d = 0; d < n_part; ++d) {
    if (d != 0) {
      // prepare the Dirichlet priors
      obtm->prepare_part();
    }
    Rcpp::Rcout << "Processing from doc " << arrl[d] <<  " to " << arru[d] << endl;
       for (int idx = arrl[d]; idx < arru[d]; idx++){
        line = Rcpp::as<std::string>(x[idx]);
        context_id = Rcpp::as<std::string>(doc_ids[idx]);
        Doc doc(line);
        if(biterms_size > 0){
          context_id_biterm = Rcpp::as<Rcpp::List>(biterms(context_id));
          term1 = context_id_biterm["term1"];
          term2 = context_id_biterm["term2"];
          cooc = context_id_biterm["cooc"];
          for (int j = 0; j < term1.size(); j++){
            if(cooc.length() > 0){
              for (int nr = 0; nr < cooc[j]; nr++){
                obtm->bs.push_back( Biterm(term1[j], term2[j]) );
              }
            }else{
              obtm->bs.push_back( Biterm(term1[j], term2[j]) );
            }
          }
        }else{
          // gen biterms
          doc.gen_biterms(obtm->bs, win);
        }
        for (int i = 0; i < doc.size(); ++i) {
          int w = doc.get_w(i);
          // the background word distribution
          obtm->pw_b[w] += 1;
        }
        // Rcpp::Rcout << "n(biterms)=" << obtm->bs.size() << endl;
      }
      obtm->pw_b.normalize();
      obtm->model_init();
      // init model
      // for (vector<Biterm>::iterator b = obtm->bs.begin(); b != obtm->bs.end(); ++b) {
      //   int k = Sampler::uni_sample(K);
      //   obtm->assign_biterm_topic(*b, k);
      // }
      double inf = std::numeric_limits<double>::infinity();
      double loglik_old = -inf;
      for (int it = 1; it < iter + 1; ++it) {
        for (unsigned int b = 0; b < obtm->bs.size(); ++b) {
           obtm->update_biterm(obtm->bs[b]);
        }
        if(trace > 0){
          if ((it) % trace == 0){
            Rcpp::Rcout << Rcpp::as<std::string>(format_posixct(sys_time())) << " End of GS iteration " << it << "/" << iter << endl;
          }
        }
        if(check_convergence > 0){
          if ((it) % check_convergence == 0){
            double loglik = obtm -> loglik();
            Rcpp::Rcout << " Loglik: " <<  loglik << endl;
            if (loglik_old/loglik - 1 < convergence_tol) {
              Rcpp::Rcout << " Achieved convergence after " << it << "/" << iter << " iterations " << endl;
              break;
            }
            loglik_old = loglik;
          }
        }
        
         Rcpp::checkUserInterrupt();
      }
      
    }
  // p(z) is determinated by the overall proportions
  // of biterms in it
  Pvec<double> pz(K);	          // p(z) = theta
  for (int k = 0; k < K; k++){
    pz[k] = (obtm->nb_z[k] + obtm->alpha[k]);
  }
  pz.normalize();

  std::vector<double> p_z;
  for (unsigned int i = 0; i < pz.size(); ++i){
    p_z.push_back(pz[i]);
  }

   Pmat<double> pw_z = obtm->nwz.toDouble() + obtm->beta;   // p(w|z) = phi, size K * M
   pw_z.normr();

  Rcpp::NumericMatrix pwz(pw_z.cols(), pw_z.rows());
    for (int m = 0; m < pw_z.rows(); ++m) {
      for (int n = 0; n < pw_z.cols(); ++n) {
        pwz(n, m) = pw_z[m][n];
      }
    }

 // Rcpp::Rcout << "n(biterms)=" << obtm->bs.size() << endl;

  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("model") = obtm,
    Rcpp::Named("K") = K,
    Rcpp::Named("W") = W,
    Rcpp::Named("alpha") = a,
    Rcpp::Named("beta") = b,
    Rcpp::Named("iter") = iter,
    Rcpp::Named("lambda") = lam,
    Rcpp::Named("theta") = p_z,
    Rcpp::Named("phi") = pwz
  );
  return out;
}

// [[Rcpp::export]]
Rcpp::List obtm_biterms(SEXP model) {
  Rcpp::XPtr<OBTM> obtm(model);
  unsigned int nr_biterms = obtm->bs.size();
  std::vector<int> term1;
  std::vector<int> term2;
  std::vector<int> topic;
  for (unsigned int i = 0; i < nr_biterms; i++){
    term1.push_back(obtm->bs[i].get_wi() + 1);
    term2.push_back(obtm->bs[i].get_wj() + 1);
    topic.push_back(obtm->bs[i].get_z() + 1);
  }
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("n") = nr_biterms,
    Rcpp::Named("biterms") = Rcpp::List::create(
      Rcpp::Named("term1") = term1,
      Rcpp::Named("term2") = term2,
      Rcpp::Named("topic") = topic
    )
  );
  return out;
}

