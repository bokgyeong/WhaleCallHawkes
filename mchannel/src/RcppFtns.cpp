// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include <RcppArmadillo.h>
#include <limits>
#include <omp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


#define minimum(x,y) (((x) < (y)) ? (x) : (y))

using namespace std;
using namespace Rcpp;
using namespace arma;



// ============================================================================-
// Compute background rate or its integral
// ============================================================================-

// Compute lam0 for HP k
// [[Rcpp::export]]
vec compLam0k(int k, vec ts, vec marks, double maxT, double beta0k, mat Xmk, vec betak, double deltak, vec Wm, vec indlam0, vec knts){
  uvec indmk = find( (marks == k) );
  int nk = indmk.size(), m = Wm.size()-1, dummy;
  vec lam0 = zeros(nk);
  vec lam0m = exp( beta0k + Xmk * betak + deltak * Wm );
  
  for(int i = 0; i < nk; i ++){
    dummy = indlam0[ indmk[i] ];
    lam0[i] = lam0m[dummy] + ( (lam0m[dummy+1] - lam0m[dummy]) / (maxT / m) ) * (ts[ indmk[i] ] - knts[dummy]);
  }
  
  return lam0;
}



// [[Rcpp::export]]
vec compLam0k_parallel(int k, vec ts, vec marks, double maxT, double beta0k, mat Xmk, vec betak, double deltak, vec Wm, vec indlam0, vec knts, int ncores){
  uvec indmk = find( (marks == k) );
  int nk = indmk.size(), m = Wm.size()-1;
  vec lam0 = zeros(nk);
  vec lam0m = exp( beta0k + Xmk * betak + deltak * Wm );
  
#pragma omp parallel num_threads(ncores)
{
#pragma omp for
  for(int i = 0; i < nk; i ++){
    int dummy = indlam0[ indmk[i] ];
    lam0[i] = lam0m[dummy] + ( (lam0m[dummy+1] - lam0m[dummy]) / (maxT / m) ) * (ts[ indmk[i] ] - knts[dummy]);
  }
}

return lam0;
}



// Compute lam0
// [[Rcpp::export]]
vec compLam0(vec ts, vec marks, double maxT, vec beta0, List Xm, mat beta, vec delta, vec Wm, vec indlam0, vec knts){
  int n = ts.size(), m = Wm.size()-1, K = beta0.size(), dummy;
  vec lam0m(m+1), lam0 = zeros(n);
  
  for(int i = 0; i < n; i ++){
    dummy = indlam0[i];
    lam0m = exp( beta0[ marks[i] ] + as<mat>(Xm[marks[i]]) * beta.col( marks[i] ) + delta[ marks[i] ] * Wm );
    lam0[i] = lam0m[dummy] + ( (lam0m[dummy+1] - lam0m[dummy]) / (maxT / m) ) * (ts[i] - knts[dummy]);
  }
  
  return lam0;
}



// [[Rcpp::export]]
vec compLam0_parallel(vec ts, vec marks, double maxT, vec beta0, List Xm, mat beta, vec delta, vec Wm, vec indlam0, vec knts, int ncores){
  int n = ts.size(), m = Wm.size()-1, K = beta0.size();
  vec lam0 = zeros(n);
  
#pragma omp parallel num_threads(ncores)
{
#pragma omp for
  for(int i = 0; i < n; i ++){
    int dummy = indlam0[i];
    vec lam0m = exp( beta0[ marks[i] ] + as<mat>(Xm[marks[i]]) * beta.col( marks[i] ) + delta[ marks[i] ] * Wm );
    lam0[i] = lam0m[dummy] + ( (lam0m[dummy+1] - lam0m[dummy]) / (maxT / m) ) * (ts[i] - knts[dummy]);
  }
}

return lam0;
}



// Compute lam0 at MARU k
// [[Rcpp::export]]
vec compLam0z0k(int k, vec ts, vec marks, vec branching, double maxT, double beta0k, mat Xmk, vec betak, double deltak, vec Wm, vec indlam0, vec knts){
  uvec indmk = find( (marks == k) && (branching == 0) );
  int nk = indmk.size(), m = Wm.size()-1, dummy;
  vec lam0 = zeros(nk);
  vec lam0m = exp( beta0k + Xmk * betak + deltak * Wm );
  
  for(int i = 0; i < nk; i ++){
    dummy = indlam0[ indmk[i] ];
    lam0[i] = lam0m[dummy] + ( (lam0m[dummy+1] - lam0m[dummy]) / (maxT / m) ) * (ts[ indmk[i] ] - knts[dummy]);
  }
  
  return lam0;
}



// [[Rcpp::export]]
vec compLam0z0k_parallel(int k, vec ts, vec marks, vec branching, double maxT, double beta0k, mat Xmk, vec betak, double deltak, vec Wm, vec indlam0, vec knts, int ncores){
  uvec indmk = find( (marks == k) && (branching == 0) );
  int nk = indmk.size(), m = Wm.size()-1;
  vec lam0 = zeros(nk);
  vec lam0m = exp( beta0k + Xmk * betak + deltak * Wm );
  
#pragma omp parallel num_threads(ncores)
{
#pragma omp for
  for(int i = 0; i < nk; i ++){
    int dummy = indlam0[ indmk[i] ];
    lam0[i] = lam0m[dummy] + ( (lam0m[dummy+1] - lam0m[dummy]) / (maxT / m) ) * (ts[ indmk[i] ] - knts[dummy]);
  }
}

return lam0;
}



// Approximation to integral of lam0 for MARU k
// [[Rcpp::export]]
double compIntLam0k(double maxT, double beta0k, mat Xmk, vec betak, double deltak, vec Wm){
  int m = Wm.size()-1;
  vec intLam0 = zeros(m);
  vec lam0m = exp( beta0k + Xmk * betak + deltak * Wm );
  
  for(int i = 0; i < m; i ++){
    intLam0[i] = (lam0m[i] + lam0m[i+1]) * maxT / (2 * m);
  }
  
  return sum(intLam0);
}



// [[Rcpp::export]]
double compIntLam0k_parallel(double maxT, double beta0k, mat Xmk, vec betak, double deltak, vec Wm, int ncores){
  int m = Wm.size()-1;
  vec intLam0 = zeros(m);
  vec lam0m = exp( beta0k + Xmk * betak + deltak * Wm );
  
#pragma omp parallel num_threads(ncores)
{
#pragma omp for
  for(int i = 0; i < m; i ++){
    intLam0[i] = (lam0m[i] + lam0m[i+1]) * maxT / (2 * m);
  }
}
return sum(intLam0);
}



// Approximation to integral of lam0
// [[Rcpp::export]]
double compIntLam0(double maxT, vec beta0, List Xm, mat beta, vec delta, vec Wm){
  int m = Wm.size()-1, K = beta0.size();
  vec lam0m(m+1);
  mat intLam0 = zeros(m, K);
  
  for(int k = 0; k < K; k ++){
    lam0m = exp( beta0[k] + as<mat>(Xm[k]) * beta.col(k) + delta[k] * Wm );
    for(int i = 0; i < m; i ++){
      intLam0(i,k) = (lam0m[i] + lam0m[i+1]) * maxT / (2 * m);
    }
  }
  
  return sum(sum(intLam0));
}


// [[Rcpp::export]]
double compIntLam0_parallel(double maxT, vec beta0, List Xm, mat beta, vec delta, vec Wm, int ncores){
  int m = Wm.size()-1, K = beta0.size();
  mat intLam0 = zeros(m, K);
  
#pragma omp parallel num_threads(ncores)
{
#pragma omp for
  for(int k = 0; k < K; k ++){
    vec lam0m = exp( beta0[k] + as<mat>(Xm[k]) * beta.col(k) + delta[k] * Wm );
    for(int i = 0; i < m; i ++){
      intLam0(i,k) = (lam0m[i] + lam0m[i+1]) * maxT / (2 * m);
    }
  }
}

return sum(sum(intLam0));
}



// ============================================================================-
// Compute triggering rate and its integral ----
// ============================================================================-

// Decay function
// [[Rcpp::export]]
double compH(double tdiff, double sdiff, double eta, double phi){
  return exp(-eta * tdiff) * exp(-phi * sdiff);
}


// Integral of decay function
// [[Rcpp::export]]
double compIntH(double tdiff, double sdiff, double eta, double phi){
  return  exp(-phi * sdiff) * (1 - exp(-eta * tdiff)) / eta;
}


// [[Rcpp::export]]
double compSumIntHl(int l, vec ts, vec marks, double maxT, mat distmat, double eta, double phi){
  uvec dummy = find( (marks == l) );
  vec tsl = ts.elem( dummy );
  int nl = tsl.size(), K = distmat.n_rows;
  mat dummyrate = zeros(nl, K);
  
  for(int i = 0; i < nl; i ++){
    for(int k = 0; k < K; k ++){
      dummyrate(i,k) = compIntH(maxT - tsl[i], distmat(l, k), eta, phi);
    }
  }
  return sum(sum(dummyrate));
}



// [[Rcpp::export]]
double compSumIntHl_parallel(int l, vec ts, vec marks, double maxT, mat distmat, double eta, double phi, int ncores){
  uvec dummy = find( (marks == l) );
  vec tsl = ts.elem( dummy );
  int nl = tsl.size(), K = distmat.n_rows;
  mat dummyrate = zeros(nl, K);
  
#pragma omp parallel num_threads(ncores)
{
#pragma omp for
  for(int i = 0; i < nl; i ++){
    for(int k = 0; k < K; k ++){
      dummyrate(i,k) = compIntH(maxT - tsl[i], distmat(l, k), eta, phi);
    }
  }
}

return sum(sum(dummyrate));
}



// [[Rcpp::export]]
double compSumIntH(vec ts, vec marks, double maxT, mat distmat, vec zeta, double eta, double phi){
  int n = ts.size(), K = distmat.n_cols;
  mat dummyrate = zeros(n, K);
  
  for(int i = 0; i < n; i ++){
    for(int k = 0; k < K; k ++){
      dummyrate(i,k) = zeta[ marks[i] ] * compIntH(maxT - ts[i], distmat(marks[i], k), eta, phi);
    }
  }
  
  return sum(sum(dummyrate));
}



// [[Rcpp::export]]
double compSumIntH_parallel(vec ts, vec marks, double maxT, mat distmat, vec zeta, double eta, double phi, int ncores){
  int n = ts.size(), K = distmat.n_cols;
  mat dummyrate = zeros(n, K);
  
#pragma omp parallel num_threads(ncores)
{
#pragma omp for
  for(int i = 0; i < n; i ++){
    for(int k = 0; k < K; k ++){
      dummyrate(i,k) = zeta[ marks[i] ] * compIntH(maxT - ts[i], distmat(marks[i], k), eta, phi);
    }
  }
}

return sum(sum(dummyrate));
}



// ============================================================================-
// Random time change theorem (RTCT) ----
// ============================================================================-

// [[Rcpp::export]]
double rtctIntLam0(double tsi, int indlam0i, vec lam0m, double maxT, vec knts){
  vec intlam0i = zeros(indlam0i+1);
  int m = lam0m.size()-1;
  
  for(int j = 0; j < indlam0i; j ++){
    intlam0i[j] = (lam0m[j] + lam0m[j+1]) * maxT / (2 * m);
  }
  
  double lam0i = lam0m[ indlam0i ] +
    ( (lam0m[ indlam0i + 1 ] - lam0m[ indlam0i ]) / (maxT / m) ) *
    (tsi - knts[ indlam0i ]);
  
  intlam0i[indlam0i] = (lam0m[ indlam0i ] + lam0i) * (tsi - knts[ indlam0i ]) / 2;
  
  return sum(intlam0i);
}



// [[Rcpp::export]]
mat rtctIntLam0ik(double tsi, int indlam0i, double maxT, vec knts, List Xm, mat beta0, List beta, mat delta, mat Wm){
  int niters = beta0.n_rows, K = beta0.n_cols;
  vec lam0m;
  mat intlam0i = zeros(niters, K);
  
  for(int k = 0; k < K; k ++){
    mat Xmk = as<mat>(Xm[k]);
    mat betak = as<mat>(beta[k]);
    for(int iter = 0; iter < niters; iter ++){
      lam0m = exp( beta0(iter,k) + Xmk * trans(betak.row(iter)) + delta(iter,k) * trans(Wm.row(iter)) );
      intlam0i(iter,k) = rtctIntLam0(tsi, indlam0i, lam0m, maxT, knts);
    }
  }
  
  return intlam0i;
}


// [[Rcpp::export]]
double rtctAlphaSumIntH(double tsi, int mi, vec ts, vec marks, mat distmat, vec alpha, double eta, double phi){
  double res = 0;
  uvec ind = find( ts < tsi );
  int sizeind = ind.size();
  
  for(int j = 0; j < sizeind; j ++){
    res = res + alpha[ marks[ ind[j] ] ] * compIntH(tsi - ts[ ind[j] ], distmat(mi, marks[ ind[j] ]), eta, phi);
  }
  
  return res;
}



// [[Rcpp::export]]
mat rtctAlphaSumIntHik(double tsi, vec ts, vec marks, mat distmat, mat alpha, vec eta, vec phi){
  int niters = eta.size(), K = distmat.n_rows;
  mat intTrigi = zeros(niters, K);
  vec alphak(K);
  double etak, phik;
  
  for(int k = 0; k < K; k ++){
    for(int iter = 0; iter < niters; iter ++){
      alphak = trans(alpha.row(iter));
      etak = eta[iter];
      phik = phi[iter];
      intTrigi(iter,k) = rtctAlphaSumIntH(tsi, k, ts, marks, distmat, alphak, etak, phik);
    }
  }
  
  return intTrigi;
}




// ============================================================================-
// Compute log-likelihood
// ============================================================================-

// Compute log-likelihood for a single posterior sample
// [[Rcpp::export]]
double compLogLiki(vec ts, vec marks, mat distmat, double maxT, mat lam0m, vec indlam0, vec knts, vec alpha, double eta, double phi){
  int m = lam0m.n_rows-1, n = ts.size(), K = distmat.n_rows, sizeind, dummy;
  double mui;
  vec lam = zeros(n), IntH = zeros(K), aSumIntH = zeros(n), lam0mk = zeros(m+1);
  mat intLam0 = zeros(m, K);
  
  for(int i = 0; i < n; i ++){
    
    // background rate
    dummy = indlam0[i];
    lam0mk = lam0m.col( marks[i] );
    mui = lam0mk[dummy] + ( (lam0mk[dummy+1] - lam0mk[dummy]) / (maxT / m) ) * (ts[i] - knts[dummy]);
    
    if( sum(alpha) == 0 ){
      // conditional intensity
      lam[i] = mui;
      
    } else {
      // self-exciting rate
      uvec ind = find( ts < ts[i] );
      sizeind = ind.size();
      vec ahi = zeros(sizeind);
      for(int j = 0; j < sizeind; j ++){
        ahi[j] = alpha[marks[j]] * exp( -eta * (ts[i] - ts[j])) * exp( -phi * distmat(marks[i], marks[j]));
      }
      
      // conditional intensity
      lam[i] = mui + sum(ahi);
      
      // integral of self-exciting rate
      for(int k = 0; k < K; k ++){
        IntH[k] = compIntH(maxT - ts[i], distmat(k, marks[i]), eta, phi);
      }
      aSumIntH[i] = alpha[marks[i]] * sum(IntH);
    }
  }
  
  // integral of background rate
  for(int k = 0; k < K; k ++){
    lam0mk = lam0m.col( k );
    for(int i = 0; i < m; i ++){
      intLam0(i,k) = (lam0mk[i] + lam0mk[i+1]) * maxT / (2 * m);
    }
  }
  
  // loglikelihood
  if( sum(alpha) == 0 ){
    return sum( log(lam) ) - sum(sum(intLam0));
    
  } else {
    return sum( log(lam) ) - sum(sum(intLam0)) - sum(aSumIntH);
  }
}



// Compute posterior mean and CI of the intensity at time t 
// [[Rcpp::export]]
List compLamt(double t, vec ts, vec marks, mat distmat, vec knts, double maxT, double indlam0t, 
              List Xm, mat beta0, List beta, mat delta, mat Wm, mat alpha, vec eta, vec phi){
  
  int niters = beta0.n_rows, m = Wm.n_cols, K = distmat.n_rows, sizeind;
  vec lam0m;
  mat Xmk, betak;
  mat muk = zeros(niters, K), sek = zeros(niters, K);
  vec P = { 0.025, 0.5, 0.975 };
  mat resBack = zeros(3, 11), resSE = zeros(3, 11), resLam = zeros(3, 11);
  vec alphaiter(K);
  double etaiter, phiiter;
  
  // background rate
  for(int k = 0; k < K; k ++){
    Xmk = as<mat>(Xm[k]);
    betak = as<mat>(beta[k]);
    
    for(int iter = 0; iter < niters; iter ++){
      lam0m = exp( beta0(iter,k) + Xmk * trans(betak.row(iter)) + delta(iter,k) * trans(Wm.row(iter)) );
      muk(iter, k) = lam0m[indlam0t] + ( (lam0m[indlam0t+1] - lam0m[indlam0t]) / (maxT / m) ) * (t - knts[indlam0t]);  
    }
    resBack.col(k) = quantile(muk.col(k), P);
  }
  resBack.col(K) = quantile(sum(muk, 1), P);
  Rprintf("Background intensity is done\n");
  
  
  // self-exciting rate
  for(int k = 0; k < K; k ++){
    for(int iter = 0; iter < niters; iter ++){
      alphaiter = trans(alpha.row(iter));
      etaiter = eta[iter];
      phiiter = phi[iter];
      
      uvec ind = find( ts < t );
      sizeind = ind.size();
      vec ahi = zeros(sizeind);
      for(int j = 0; j < sizeind; j ++){
        ahi[j] = alphaiter[ marks[j] ] * exp( -etaiter * (t - ts[j])) * exp( -phiiter * distmat(k, marks[j]));
      }
      sek(iter,k) = sum(ahi);  
    }
    resSE.col(k) = quantile(sek.col(k), P);
  }
  resSE.col(K) = quantile(sum(sek, 1), P);
  Rprintf("Countercall intensity is done\n");
  
  
  // intensity
  for(int k = 0; k < K; k ++){
    resLam.col(k) = quantile(muk.col(k) + sek.col(k), P);
  }
  resLam.col(K) = quantile(sum(muk, 1) + sum(sek, 1), P);
  Rprintf("Intensity is done\n");
  
  
  return Rcpp::List::create(Rcpp::Named("resBack") = resBack,
                            Rcpp::Named("resSE") = resSE,
                            Rcpp::Named("resLam") = resLam);
}




// ============================================================================-
// Sample latent branching structure ----
// ============================================================================-

// [[Rcpp::export]]
vec sampleBranching(vec ts, vec marks, mat distmat, vec lam0, vec zeta, double eta, double phi) {
  int n = ts.size(), l, k, parent, ndummy;
  double temp;
  vec branching = zeros(n), probs;
  uvec dummy;
  
  std::random_device rd;
  std::mt19937 gen(rd());
  
  for (int i = 1; i < n; i++) {
    
    k = marks[i];
    dummy = find( ts < ts[i] );
    ndummy = dummy.size();
    
    if(ndummy == 0){
      branching[i] = 0;
      
    } else {
      probs = zeros(ndummy + 1);
      probs[0] = lam0[i];
      
      for (int j = 0; j < ndummy; j++) {
        l = marks[ dummy[j] ];
        temp = zeta[l] * compH(ts[i] - ts[ dummy[j] ], distmat(l,k), eta, phi);
        probs[j+1] = temp;
      }
      probs = probs / sum(probs);
      std::discrete_distribution<> d(probs.begin(), probs.end());
      
      parent = d(gen);
      branching[i] = parent;
    }
  }
  
  return branching;
}




// ============================================================================-
// Distribution functions
// ============================================================================-

// Log of unnormalized Gaussian pdf
// [[Rcpp::export]]
double normal_logh(double y, double mu, double sig2){
  return - 0.5 * pow(y - mu, 2) / sig2;
}


// Log of unnormalized Gamma pdf
// [[Rcpp::export]]
double gamma_logh(double y, double alpha, double beta){
  return (alpha-1) * log(y) - beta * y;
}


// Log of unnormalized inverse Gamma pdf
// [[Rcpp::export]]
double invgamma_logh(double y, double alpha, double beta){
  return (-alpha-1) * log(y) - beta / y;
}


// Logo of unnormalized MVN pdf
// [[Rcpp::export]]
double MVN_logh(vec y, vec mu, mat invSigma){
  vec result = - 0.5 * trans(y - mu) * invSigma * (y - mu);
  return result[0];
}



// ============================================================================-
// Model fitting
// ============================================================================-



// MVN prior for beta0 and beta
// invGamma for tau
// [[Rcpp::export]]
List fitNHPP(int niter, vec ts, vec marks, List Xm, double maxT,
             double rho_beta, mat distmat, vec knts,
             double betaTilde0, vec beta0, vec betaTilde, mat beta, vec tau, vec indlam0,
             double shape_tau, double rate_tau,
             vec sigma2, mat COVbeta0, List COVbeta,
             bool updateCOV, int adaptInterval, double adaptFactorExponent,
             vec adapIter, int ncores) {
  
  // double negativeInf = -std::numeric_limits<float>::infinity();;
  // double positiveInf = std::numeric_limits<float>::infinity();;
  
  int n = ts.size(), p = beta.n_rows, K = distmat.n_rows, m = knts.size() - 1;
  double logprob;
  mat postBeta0 = zeros(niter, K), postTau = zeros(niter, 1+p), postTilde = zeros(niter,p+1);
  mat postBranching = zeros(niter, n), accprob = zeros(niter, 1+p);
  List postBeta(p);
  for(int i = 0; i < p; i ++){
    postBeta[i] = zeros(niter, K);
  }
  mat postBetai = zeros(niter, K);
  double sTilde, mTilde;
  vec postj;
  vec rhat = zeros(1+p), gamma1 = zeros(1+p), gamma2 = zeros(1+p), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, shape, scale;
  vec newbeta0(K), betai(K), newbetai(K);
  mat newbeta(p, K);
  uvec dummy;
  
  mat cholCOVbeta0(K, K);
  List cholCOVbeta(p);
  cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
  for(int i = 0; i < p; i ++){
    cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
  }
  
  mat Vmat = exp(-distmat / rho_beta);
  mat invVmat = inv(Vmat), cholVmat = chol(Vmat);
  double oneinvVmatone = sum( trans(ones(K)) * invVmat * ones(K) );
  vec invVmatone = invVmat * ones(K);
  
  double intLam0, newintLam0;
  vec lam0, newlam0;
  if(ncores == 1){
    intLam0 = compIntLam0(maxT, beta0, Xm, beta, zeros(K), zeros(m+1));
    lam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts);
    
  } else {
    intLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, zeros(K), zeros(m+1), ncores);
    lam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts, ncores);
  }
  
  
  
  
  // Start MCMC
  for(int s = 0; s < niter; s++) {
    
    // M-H update for beta0
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COVbeta0 = COVbeta0 + gamma1[0] * ( cov( postBeta0.rows(s+1-adaptInterval, s-1) ) - COVbeta0 );
        cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
        
        adapIter[0] = adapIter[0] + 1;
      }
    }
    
    newbeta0 = beta0 + cholCOVbeta0 * randn(K);
    
    if(ncores == 1){
      newlam0 = compLam0(ts, marks, maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts);
      newintLam0 = compIntLam0(maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1));
      
    } else {
      newlam0 = compLam0_parallel(ts, marks, maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts, ncores);
      newintLam0 = compIntLam0_parallel(maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1), ncores);
    }
    
    
    logprob = (sum(log(newlam0)) - sum(log(lam0))) - (newintLam0 - intLam0) +
      MVN_logh(newbeta0, betaTilde0 * ones(K), invVmat/tau[0]) - MVN_logh(beta0, betaTilde0 * ones(K), invVmat/tau[0]);
    
    if (log(randu()) < logprob) {
      beta0 = newbeta0;
      lam0 = newlam0;
      intLam0 = newintLam0;
      
      accprob(s, 0) = 1;
    }
    
    postBeta0.row(s) = trans(beta0);
    if(s == 0){ Rprintf("sampled beta0\n"); }
    
    
    
    // Gibbs update for betaTilde0
    sTilde = 1 / (oneinvVmatone / tau[0] + 1/100);
    mTilde = sTilde / tau[0] * sum(trans(beta0) * invVmatone);
    betaTilde0 = mTilde + sqrt(sTilde) * randn();
    postTilde(s,0) = betaTilde0;
    if(s == 0){ Rprintf("sampled betaTilde0\n"); }
    
    
    
    // M-H updates for beta
    for(int i = 0; i < p; i ++){
      if( updateCOV ){
        if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(1+i);
          rhat[1+i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
          gamma1[1+i] = 1 / pow(adapIter[1+i], c1);
          gamma2[1+i] = c0 * gamma1[1+i];
          sigma2[1+i] = exp( log(sigma2[1+i]) + gamma2[1+i] * (rhat[1+i] - ropt) );
          
          postBetai = as<mat>(postBeta[i]);
          COVbeta[i] = as<mat>(COVbeta[i]) + gamma1[1+i] * ( cov( postBetai.rows(s+1-adaptInterval, s-1) ) - as<mat>(COVbeta[i]) );
          cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
          
          adapIter[1+i] = adapIter[1+i] + 1;
        }
      }
      
      betai = trans(beta.row(i));
      newbetai = betai + as<mat>(cholCOVbeta[i]) * randn(K);
      newbeta = beta;
      newbeta.row(i) = trans(newbetai);
      
      if(ncores == 1){
        newlam0 = compLam0(ts, marks, maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1), indlam0, knts);
        newintLam0 = compIntLam0(maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1));
        
      } else {
        newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1), indlam0, knts, ncores);
        newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1), ncores);
      }
      
      logprob = (sum(log(newlam0)) - sum(log(lam0))) - (newintLam0 - intLam0) +
        MVN_logh(newbetai, betaTilde[i] * zeros(K), invVmat / tau[1+i]) - MVN_logh(betai, betaTilde[i] * zeros(K), invVmat / tau[1+i]);
      
      if (log(randu()) < logprob) {
        betai = newbetai;
        beta = newbeta;
        lam0 = newlam0;
        intLam0 = newintLam0;
        
        accprob(s, 1+i) = 1;
      }
      
      postBetai = as<mat>(postBeta[i]);
      postBetai.row(s) = trans(betai);
      postBeta[i] = postBetai;
    }
    if(s == 0){ Rprintf("sampled beta\n"); }
    
    
    
    // Gibbs update for betaTilde
    for(int i = 0; i < p; i ++){
      sTilde = 1 / (oneinvVmatone / tau[1+i] + 1/100);
      mTilde = sTilde / tau[1+i] * sum(beta.row(i) * invVmatone);
      betaTilde[i] = mTilde + sqrt(sTilde) * randn();
      postTilde(s,1+i) = betaTilde[i]; 
    }
    if(s == 0){ Rprintf("sampled betaTilde\n"); }
    
    
    
    // Gibbs update for tau
    shape = shape_tau + 0.5 * K;
    scale = ( 1 / ( rate_tau + 0.5 * trans(beta0 - betaTilde0) * invVmat * (beta0 - betaTilde0) ) )[0];
    tau[0] = 1/randg(distr_param(shape, scale));
    
    for(int i = 0; i < p; i ++){
      shape = shape_tau + 0.5 * K;
      scale = ( 1 / ( rate_tau + 0.5 * (beta.row(i) - betaTilde[i]) * invVmat * trans(beta.row(i) - betaTilde[i]) ) )[0];
      
      tau[1+i] = 1/randg(distr_param(shape, scale));
    }
    
    postTau.row(s) = trans(tau);
    if(s == 0){ Rprintf("sampled tau\n"); }
    
    
    
    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postBeta0") = postBeta0,
                            Rcpp::Named("postBeta") = postBeta,
                            Rcpp::Named("postTilde") = postTilde,
                            Rcpp::Named("postTau") = postTau,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COVbeta0") = COVbeta0,
                            Rcpp::Named("COVbeta") = COVbeta);
}




// ESS for W, fixed rho, variance = 1, mean = 0
// MVN prior for beta0, beta, and delta
// invGamma for tau
// [[Rcpp::export]]
List fitLGCP(int niter, vec ts, vec marks, List Xm, double maxT,
             double rho_beta, double rho_w,
             mat distmat, vec knts, mat tdiffm,
             double betaTilde0, vec beta0, vec betaTilde, mat beta, vec tau,
             vec Wm, double deltaTilde, vec delta, vec indlam0,
             double shape_tau, double rate_tau,
             vec sigma2, mat COVbeta0, List COVbeta, mat COVdelta,
             bool updateCOV, int adaptInterval, double adaptFactorExponent,
             vec adapIter, int ncores) {
  
  double negativeInf = -std::numeric_limits<float>::infinity();;
  // double positiveInf = std::numeric_limits<float>::infinity();;
  
  int n = ts.size(), p = beta.n_rows, K = distmat.n_rows, m = knts.size() - 1;
  double logprob;
  mat postBeta0 = zeros(niter, K), postWm(niter, (m+1)), postDelta(niter, K), postTau = zeros(niter, 1+p+1), postTilde = zeros(niter,p+2);
  mat accprob = zeros(niter, 1+p+1);
  List postBeta(p);
  for(int i = 0; i < p; i ++){
    postBeta[i] = zeros(niter, K);
  }
  mat postBetai = zeros(niter, K);
  vec postj;
  double sTilde, mTilde;
  int count;
  vec newWm(m+1);
  vec rhat = zeros(1+p+1), gamma1 = zeros(1+p+1), gamma2 = zeros(1+p+1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, shape, scale;
  vec xvals;
  double llprev, llthred, llnew, thetamin, thetamax, theta;
  vec priorWm, priorbeta0;
  bool accept;
  vec newbeta0(K), betai(K), newbetai(K);
  mat newbeta(p, K);
  vec newdelta(K);
  
  mat cholCOVbeta0(K, K), cholCOVdelta(K, K);
  List cholCOVbeta(p);
  cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
  for(int i = 0; i < p; i ++){
    cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
  }
  cholCOVdelta = trans( chol( sigma2[1+p] * ( COVdelta + 0.0000000001 * diagmat(ones(K)) ) ) );
  
  mat Vmat = exp(-distmat / rho_beta);
  mat Sigma = exp(-tdiffm / rho_w);
  mat invSigma = inv(Sigma), invVmat = inv(Vmat), cholSigma = chol(Sigma), cholVmat = chol(Vmat);
  double oneinvVmatone = sum( trans(ones(K)) * invVmat * ones(K) );
  vec invVmatone = invVmat * ones(K);
  
  double intLam0, newintLam0;
  vec lam0, newlam0;
  if(ncores == 1){
    intLam0 = compIntLam0(maxT, beta0, Xm, beta, exp(delta), Wm);
    lam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, exp(delta), Wm, indlam0, knts);
    
  } else {
    intLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, exp(delta), Wm, ncores);
    lam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, exp(delta), Wm, indlam0, knts, ncores);
  }
  
  
  
  
  // Start MCMC
  for(int s = 0; s < niter; s++) {
    
    // M-H update for beta0
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COVbeta0 = COVbeta0 + gamma1[0] * ( cov( postBeta0.rows(s+1-adaptInterval, s-1) ) - COVbeta0 );
        cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
        
        adapIter[0] = adapIter[0] + 1;
      }
    }
    
    newbeta0 = beta0 + cholCOVbeta0 * randn(K);
    if(ncores == 1){
      newlam0 = compLam0(ts, marks, maxT, newbeta0, Xm, beta, exp(delta), Wm, indlam0, knts);
      newintLam0 = compIntLam0(maxT, newbeta0, Xm, beta, exp(delta), Wm);
      
    } else {
      newlam0 = compLam0_parallel(ts, marks, maxT, newbeta0, Xm, beta, exp(delta), Wm, indlam0, knts, ncores);
      newintLam0 = compIntLam0_parallel(maxT, newbeta0, Xm, beta, exp(delta), Wm, ncores);
    }
    
    logprob = (sum(log(newlam0)) - sum(log(lam0))) - (newintLam0 - intLam0) +
      MVN_logh(newbeta0, betaTilde0 * ones(K), invVmat/tau[0]) - MVN_logh(beta0, betaTilde0 * ones(K), invVmat/tau[0]);
    
    if (log(randu()) < logprob) {
      beta0 = newbeta0;
      lam0 = newlam0;
      intLam0 = newintLam0;
      lam0 = newlam0;
      
      accprob(s, 0) = 1;
    }
    
    postBeta0.row(s) = trans(beta0);
    if(s == 0){ Rprintf("sampled beta0\n"); }
    
    
    
    // Gibbs update for betaTilde0
    sTilde = 1 / (oneinvVmatone / tau[0] + 1/100);
    mTilde = sTilde / tau[0] * sum(trans(beta0) * invVmatone);
    betaTilde0 = mTilde + sqrt(sTilde) * randn();
    postTilde(s,0) = betaTilde0;
    if(s == 0){ Rprintf("sampled betaTilde0\n"); }
    
    
    
    // M-H updates for beta
    for(int i = 0; i < p; i ++){
      if( updateCOV ){
        if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(1+i);
          rhat[1+i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
          gamma1[1+i] = 1 / pow(adapIter[1+i], c1);
          gamma2[1+i] = c0 * gamma1[1+i];
          sigma2[1+i] = exp( log(sigma2[1+i]) + gamma2[1+i] * (rhat[1+i] - ropt) );
          
          postBetai = as<mat>(postBeta[i]);
          COVbeta[i] = as<mat>(COVbeta[i]) + gamma1[1+i] * ( cov( postBetai.rows(s+1-adaptInterval, s-1) ) - as<mat>(COVbeta[i]) );
          cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
          
          adapIter[1+i] = adapIter[1+i] + 1;
        }
      }
      
      betai = trans(beta.row(i));
      newbetai = betai + as<mat>(cholCOVbeta[i]) * randn(K);
      newbeta = beta;
      newbeta.row(i) = trans(newbetai);
      
      if(ncores == 1){
        newlam0 = compLam0(ts, marks, maxT, beta0, Xm, newbeta, exp(delta), Wm, indlam0, knts);
        newintLam0 = compIntLam0(maxT, beta0, Xm, newbeta, exp(delta), Wm);
        
      } else {
        newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, newbeta, exp(delta), Wm, indlam0, knts, ncores);
        newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, newbeta, exp(delta), Wm, ncores);
      }
      
      logprob = (sum(log(newlam0)) - sum(log(lam0))) - (newintLam0 - intLam0) +
        MVN_logh(newbetai, betaTilde[i] * zeros(K), invVmat / tau[1+i]) - MVN_logh(betai, betaTilde[i] * zeros(K), invVmat / tau[1+i]);
      
      if (log(randu()) < logprob) {
        betai = newbetai;
        beta = newbeta;
        lam0 = newlam0;
        intLam0 = newintLam0;
        lam0 = newlam0;
        
        accprob(s, 1+i) = 1;
      }
      
      postBetai = as<mat>(postBeta[i]);
      postBetai.row(s) = trans(betai);
      postBeta[i] = postBetai;
    }
    if(s == 0){ Rprintf("sampled beta\n"); }
    
    
    
    // Gibbs update for betaTilde
    for(int i = 0; i < p; i ++){
      sTilde = 1 / (oneinvVmatone / tau[1+i] + 1/100);
      mTilde = sTilde / tau[1+i] * sum(beta.row(i) * invVmatone);
      betaTilde[i] = mTilde + sqrt(sTilde) * randn();
      postTilde(s,1+i) = betaTilde[i]; 
    }
    if(s == 0){ Rprintf("sampled betaTilde\n"); }
    
    
    
    // Elliptical slice sampling for Wm
    llprev = sum(log(lam0)) - intLam0;
    priorWm = trans( cholSigma ) * randn(m+1);
    
    thetamin = 0;
    thetamax = 2 * M_PI;
    theta = thetamin + randu() * (thetamax - thetamin);
    thetamin = theta - 2 * M_PI;
    thetamax = theta;
    
    llthred = llprev + log(randu());
    accept = false;
    count = 0;
    
    while(accept == false){
      count = count + 1;
      
      newWm = Wm * cos(theta) + priorWm * sin(theta);
      
      if(ncores == 1){
        newlam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, exp(delta), newWm, indlam0, knts);
        newintLam0 = compIntLam0(maxT, beta0, Xm, beta, exp(delta), newWm);
        
      } else {
        newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, exp(delta), newWm, indlam0, knts, ncores);
        newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, exp(delta), newWm, ncores);
      }
      
      llnew = sum(log(newlam0)) - newintLam0;
      
      if(llnew > llthred){
        llprev = llnew;
        accept = true;
      } else {
        if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
        theta = thetamin + randu() * (thetamax - thetamin);
        if ( (count) % 1000 == 0 ) {
          Rprintf("ESS for Wm: %d iterations...\n", count);
        }
      }
    }
    Wm = newWm;
    lam0 = newlam0;
    intLam0 = newintLam0;
    
    postWm.row(s) = trans(Wm);
    
    
    
    // MH updates for delta
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1+p);
        rhat[1+p] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1+p] = 1 / pow(adapIter[1+p], c1);
        gamma2[1+p] = c0 * gamma1[1+p];
        sigma2[1+p] = exp( log(sigma2[1+p]) + gamma2[1+p] * (rhat[1+p] - ropt) );
        
        COVdelta = COVdelta + gamma1[1+p] * ( cov( postDelta.rows(s+1-adaptInterval, s-1) ) - COVdelta );
        cholCOVdelta = trans( chol( sigma2[1+p] * ( COVdelta + 0.0000000001 * diagmat(ones(K)) ) ) );
        
        adapIter[1+p] = adapIter[1+p] + 1;
      }
    }
    
    newdelta = delta + cholCOVdelta * randn(K);
    if(ncores == 1){
      newlam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, exp(newdelta), Wm, indlam0, knts);
      newintLam0 = compIntLam0(maxT, beta0, Xm, beta, exp(newdelta), Wm);
      
    } else {
      newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, exp(newdelta), Wm, indlam0, knts, ncores);
      newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, exp(newdelta), Wm, ncores);
    }
    
    logprob = (sum(log(newlam0)) - sum(log(lam0))) - (newintLam0 - intLam0) +
      MVN_logh(newdelta, deltaTilde * ones(K), invVmat/tau[1+p]) - MVN_logh(delta, deltaTilde * ones(K), invVmat/tau[1+p]);
    
    if (log(randu()) < logprob) {
      delta = newdelta;
      lam0 = newlam0;
      intLam0 = newintLam0;
      
      accprob(s, 1+p) = 1;
    }
    
    postDelta.row(s) = trans(delta);
    if(s == 0){ Rprintf("sampled delta\n"); }
    
    
    
    // Gibbs update for deltaTilde
    sTilde = 1 / (oneinvVmatone / tau[1+p] + 1/100);
    mTilde = sTilde / tau[1+p] * sum(trans(delta) * invVmatone);
    deltaTilde = mTilde + sqrt(sTilde) * randn();
    postTilde(s,1+p) = deltaTilde;
    if(s == 0){ Rprintf("sampled deltaTilde\n"); }
    
    
    
    // Gibbs update for tau
    shape = shape_tau + 0.5 * K;
    scale = ( 1 / ( rate_tau + 0.5 * trans(beta0 - betaTilde0) * invVmat * (beta0 - betaTilde0) ) )[0];
    tau[0] = 1/randg(distr_param(shape, scale));
    
    for(int i = 0; i < p; i ++){
      shape = shape_tau + 0.5 * K;
      scale = ( 1 / ( rate_tau + 0.5 * (beta.row(i) - betaTilde[i]) * invVmat * trans(beta.row(i) - betaTilde[i]) ) )[0];
      tau[1+i] = 1/randg(distr_param(shape, scale));
    }
    
    shape = shape_tau + 0.5 * K;
    scale = ( 1 / ( rate_tau + 0.5 * trans(delta - deltaTilde) * invVmat * (delta - deltaTilde) ) )[0];
    tau[1+p] = 1/randg(distr_param(shape, scale));
    
    postTau.row(s) = trans(tau);
    if(s == 0){ Rprintf("sampled tau\n"); }
    
    
    
    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postBeta0") = postBeta0,
                            Rcpp::Named("postBeta") = postBeta,
                            Rcpp::Named("postTau") = postTau,
                            Rcpp::Named("postDelta") = postDelta,
                            Rcpp::Named("postTilde") = postTilde,
                            Rcpp::Named("postWm") = postWm,
                            Rcpp::Named("Wm") = Wm,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COVbeta0") = COVbeta0,
                            Rcpp::Named("COVbeta") = COVbeta,
                            Rcpp::Named("COVdelta") = COVdelta);
}






// MVN prior for beta0 and beta
// Gamma prior for zeta
// Unif prior for eta and phi
// [[Rcpp::export]]
List fitNHPPSE(int niter, vec ts, vec marks, List Xm, double maxT,
               double rho_beta,
               mat distmat, vec knts,
               double betaTilde0, vec beta0, vec betaTilde, mat beta, vec tau, 
               vec indlam0,
               vec zeta, double eta, double phi,
               double shape_tau, double rate_tau,
               double shape_zeta, double rate_zeta,
               double lb_eta, double ub_eta,
               double lb_phi, double ub_phi,
               vec sigma2, mat COVbeta0, List COVbeta,
               double COVeta, double COVphi,
               bool updateCOV, int adaptInterval, double adaptFactorExponent,
               vec adapIter, int ncores) {
  
  double negativeInf = -std::numeric_limits<float>::infinity();;
  // double positiveInf = std::numeric_limits<float>::infinity();;
  
  int n = ts.size(), p = beta.n_rows, K = distmat.n_rows, m = knts.size() - 1;
  double logprob;
  mat postBeta0 = zeros(niter, K), postZeta = zeros(niter, K), postTau = zeros(niter, 1+p), postTilde = zeros(niter,p+1);
  vec postEta = zeros(niter), postPhi = zeros(niter);
  mat postBranching = zeros(niter, n), accprob = zeros(niter, 1+p+1+1);
  List postBeta(p);
  for(int i = 0; i < p; i ++){
    postBeta[i] = zeros(niter, K);
  }
  mat postBetai = zeros(niter, K);
  vec postj;
  double sTilde, mTilde;
  vec branching(n), numtrig(K);
  int count, numtrigtotal, numtrigl;
  uvec indtrig, indparent, indparentl;
  double tdifftrig, sdifftrig, dummySum, dummynewSum;
  int numback;
  vec lam0z0, newlam0z0;
  vec rhat = zeros(1+p+1+1), gamma1 = zeros(1+p+1+1), gamma2 = zeros(1+p+1+1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, shape, scale;
  vec xvals;
  double llprev, llthred, llnew, thetamin, thetamax, theta;
  vec priorbeta0;
  bool accept;
  double zetal, neweta, newphi;
  vec newbeta0(K), betai(K), newbetai(K);
  mat newbeta(p, K);
  uvec indback, dummy;
  
  mat cholCOVbeta0(K, K);
  List cholCOVbeta(p);
  double cholCOVeta, cholCOVphi;
  cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
  for(int i = 0; i < p; i ++){
    cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
  }
  cholCOVeta = sqrt( sigma2[1+p] * (COVeta + 0.0000000001) );
  cholCOVphi = sqrt( sigma2[1+p+1] * (COVphi + 0.0000000001) );
  
  mat Vmat = exp(-distmat / rho_beta);
  mat invVmat = inv(Vmat), cholVmat = chol(Vmat);
  double oneinvVmatone = sum( trans(ones(K)) * invVmat * ones(K) );
  vec invVmatone = invVmat * ones(K);
  
  double intLam0, newintLam0;
  vec lam0, newlam0;
  if(ncores == 1){
    intLam0 = compIntLam0(maxT, beta0, Xm, beta, zeros(K), zeros(m+1));
    lam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts);
    
  } else {
    intLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, zeros(K), zeros(m+1), ncores);
    lam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts, ncores);
  }
  
  
  
  
  // Start MCMC
  for(int s = 0; s < niter; s++) {
    
    // Gibbs update for latent branching structure
    branching = sampleBranching(ts, marks, distmat, lam0, zeta, eta, phi);
    postBranching.row(s) = trans( branching );
    if(s == 0){ Rprintf("sampled braching structure\n"); }
    
    
    // Self-exciting process
    indtrig = find( branching > 0 ); // triggered calls
    numtrigtotal = indtrig.size();
    indparent = indtrig;
    
    tdifftrig = 0;
    sdifftrig = 0;
    
    if(numtrigtotal != 0){
      for(int i = 0; i < numtrigtotal; i ++){
        indparent[i] = branching[ indtrig[i] ] - 1; // the delayed calls' parents
      }
      
      for(int l = 0; l < K; l ++){
        indparentl = find( (marks.elem(indparent) == l) );
        numtrigl = indparentl.size();
        
        if(numtrigl == 0){
          numtrig[l] = 0;
          
        } else {
          for(int i = 0; i < numtrigl; i ++){
            tdifftrig = tdifftrig + ts[ indtrig[ indparentl[i] ] ] - ts[ indparent[ indparentl[i] ] ];
            sdifftrig = sdifftrig + distmat(l, marks[ indtrig[ indparentl[i] ] ]);
          }
          numtrig[l] = numtrigl;
        }
      }
    }
    
    
    // Background process
    indback = find(branching == 0);
    lam0z0 = lam0.elem(indback);
    numback = indback.size();
    
    
    
    // M-H update for beta0
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COVbeta0 = COVbeta0 + gamma1[0] * ( cov( postBeta0.rows(s+1-adaptInterval, s-1) ) - COVbeta0 );
        cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
        
        adapIter[0] = adapIter[0] + 1;
      }
    }
    
    newbeta0 = beta0 + cholCOVbeta0 * randn(K);
    if(ncores == 1){
      newlam0 = compLam0(ts, marks, maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts);
      newintLam0 = compIntLam0(maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1));
      
    } else {
      newlam0 = compLam0_parallel(ts, marks, maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1), indlam0, knts, ncores);
      newintLam0 = compIntLam0_parallel(maxT, newbeta0, Xm, beta, zeros(K), zeros(m+1), ncores);
    }
    newlam0z0 = newlam0.elem(indback);
    
    logprob = (sum(log(newlam0z0)) - sum(log(lam0z0))) - (newintLam0 - intLam0) +
      MVN_logh(newbeta0, betaTilde0 * ones(K), invVmat/tau[0]) - MVN_logh(beta0, betaTilde0 * ones(K), invVmat/tau[0]);
    
    if (log(randu()) < logprob) {
      beta0 = newbeta0;
      lam0 = newlam0;
      intLam0 = newintLam0;
      lam0z0 = newlam0z0;
      
      accprob(s, 0) = 1;
    }
    
    postBeta0.row(s) = trans(beta0);
    if(s == 0){ Rprintf("sampled beta0\n"); }
    
    
    
    // Gibbs update for betaTilde0
    sTilde = 1 / (oneinvVmatone / tau[0] + 1/100);
    mTilde = sTilde / tau[0] * sum(trans(beta0) * invVmatone);
    betaTilde0 = mTilde + sqrt(sTilde) * randn();
    postTilde(s,0) = betaTilde0;
    if(s == 0){ Rprintf("sampled betaTilde0\n"); }
    
    
    
    // M-H updates for beta
    for(int i = 0; i < p; i ++){
      if( updateCOV ){
        if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(1+i);
          rhat[1+i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
          gamma1[1+i] = 1 / pow(adapIter[1+i], c1);
          gamma2[1+i] = c0 * gamma1[1+i];
          sigma2[1+i] = exp( log(sigma2[1+i]) + gamma2[1+i] * (rhat[1+i] - ropt) );
          
          postBetai = as<mat>(postBeta[i]);
          COVbeta[i] = as<mat>(COVbeta[i]) + gamma1[1+i] * ( cov( postBetai.rows(s+1-adaptInterval, s-1) ) - as<mat>(COVbeta[i]) );
          cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
          
          adapIter[1+i] = adapIter[1+i] + 1;
        }
      }
      
      betai = trans(beta.row(i));
      newbetai = betai + as<mat>(cholCOVbeta[i]) * randn(K);
      newbeta = beta;
      newbeta.row(i) = trans(newbetai);
      
      if(ncores == 1){
        newlam0 = compLam0(ts, marks, maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1), indlam0, knts);
        newintLam0 = compIntLam0(maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1));
        
      } else {
        newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1), indlam0, knts, ncores);
        newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, newbeta, zeros(K), zeros(m+1), ncores);
      }
      newlam0z0 = newlam0.elem(indback);
      
      logprob = (sum(log(newlam0z0)) - sum(log(lam0z0))) - (newintLam0 - intLam0) +
        MVN_logh(newbetai, betaTilde[i] * zeros(K), invVmat / tau[1+i]) - MVN_logh(betai, betaTilde[i] * zeros(K), invVmat / tau[1+i]);
      
      if (log(randu()) < logprob) {
        betai = newbetai;
        beta = newbeta;
        lam0 = newlam0;
        intLam0 = newintLam0;
        lam0z0 = newlam0z0;
        
        accprob(s, 1+i) = 1;
      }
      
      postBetai = as<mat>(postBeta[i]);
      postBetai.row(s) = trans(betai);
      postBeta[i] = postBetai;
      
      betaTilde[i] = mean(betai); // hierarchical centering
    }
    if(s == 0){ Rprintf("sampled beta\n"); }
    
    
    
    // Gibbs update for betaTilde
    for(int i = 0; i < p; i ++){
      sTilde = 1 / (oneinvVmatone / tau[1+i] + 1/100);
      mTilde = sTilde / tau[1+i] * sum(beta.row(i) * invVmatone);
      betaTilde[i] = mTilde + sqrt(sTilde) * randn();
      postTilde(s,1+i) = betaTilde[i]; 
    }
    if(s == 0){ Rprintf("sampled betaTilde\n"); }
    
    
    
    // Gibbs update for tau
    shape = shape_tau + 0.5 * K;
    scale = ( 1 / ( rate_tau + 0.5 * trans(beta0 - betaTilde0) * invVmat * (beta0 - betaTilde0) ) )[0];
    tau[0] = 1/randg(distr_param(shape, scale));
    
    for(int i = 0; i < p; i ++){
      shape = shape_tau + 0.5 * K;
      scale = ( 1 / ( rate_tau + 0.5 * (beta.row(i) - betaTilde[i]) * invVmat * trans(beta.row(i) - betaTilde[i]) ) )[0];
      
      tau[1+i] = 1/randg(distr_param(shape, scale));
    }
    
    postTau.row(s) = trans(tau);
    if(s == 0){ Rprintf("sampled tau\n"); }
    
    
    
    // Gibbs update for zeta
    for(int l = 0; l < K; l ++){
      
      if(ncores == 1){
        dummySum = compSumIntHl(l, ts, marks, maxT, distmat, eta, phi);
      } else{
        dummySum = compSumIntHl_parallel(l, ts, marks, maxT, distmat, eta, phi, ncores);
      }
      
      shape = shape_zeta + numtrig[l];
      scale = 1 / ( rate_zeta +  dummySum );
      
      zetal = randg(distr_param(shape, scale));
      zeta[l] = zetal;
      
      postZeta(s,l) = zetal;
    }
    if(s == 0){ Rprintf("sampled zeta\n"); }
    
    
    
    // M-H update for eta
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1+p);
        rhat[1+p] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1+p] = 1 / pow(adapIter[1+p], c1);
        gamma2[1+p] = c0 * gamma1[1+p];
        sigma2[1+p] = exp( log(sigma2[1+p]) + gamma2[1+p] * (rhat[1+p] - ropt) );
        
        COVeta = COVeta + gamma1[1+p] * ( var( postEta.rows(s+1-adaptInterval, s-1) ) - COVeta );
        cholCOVeta = sqrt( sigma2[1+p] * (COVeta + 0.0000000001) );
        
        adapIter[1+p] = adapIter[1+p] + 1;
      }
    }
    
    neweta = eta + cholCOVeta * randn();
    
    if( (neweta < lb_eta) || (neweta > ub_eta) ){
      logprob = negativeInf;
      
    } else {
      if(ncores == 1){
        dummynewSum = compSumIntH(ts, marks, maxT, distmat, zeta, neweta, phi);
        dummySum = compSumIntH(ts, marks, maxT, distmat, zeta, eta, phi);
      } else{
        dummynewSum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, neweta, phi, ncores);
        dummySum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, eta, phi, ncores);
      }
      
      if( numtrigtotal == 0 ){
        logprob = - ( dummynewSum - dummySum );
        
      } else {
        logprob = - tdifftrig * (neweta - eta) - ( dummynewSum - dummySum );
      }
    }
    
    if (log(randu()) < logprob) {
      eta = neweta;
      accprob(s, 1+p) = 1;
    }
    
    postEta[s] = eta;
    if(s == 0){ Rprintf("sampled eta\n"); }
    
    
    
    // M-H update for phi
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1+p+1);
        rhat[1+p+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1+p+1] = 1 / pow(adapIter[1+p+1], c1);
        gamma2[1+p+1] = c0 * gamma1[1+p+1];
        sigma2[1+p+1] = exp( log(sigma2[1+p+1]) + gamma2[1+p+1] * (rhat[1+p+1] - ropt) );
        
        COVphi = COVphi + gamma1[1+p+1] * ( var( postPhi.rows(s+1-adaptInterval, s-1) ) - COVphi );
        cholCOVphi = sqrt( sigma2[1+p+1] * (COVphi + 0.0000000001) );
        
        adapIter[1+p+1] = adapIter[1+p+1] + 1;
      }
    }
    
    newphi = phi + cholCOVphi * randn();
    
    if( (newphi < lb_phi) || (newphi > ub_phi) ){
      logprob = negativeInf;
      
    } else {
      if(ncores == 1){
        dummynewSum = compSumIntH(ts, marks, maxT, distmat, zeta, eta, newphi);
        dummySum = compSumIntH(ts, marks, maxT, distmat, zeta, eta, phi);
      } else{
        dummynewSum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, eta, newphi, ncores);
        dummySum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, eta, phi, ncores);
      }
      
      if( numtrigtotal == 0 ){
        logprob = - ( dummynewSum - dummySum );
        
      } else {
        logprob = - sdifftrig * (newphi - phi) - ( dummynewSum - dummySum  );
      }
    }
    
    if (log(randu()) < logprob) {
      phi = newphi;
      accprob(s, 1+p+1) = 1;
    }
    
    postPhi[s] = phi;
    if(s == 0){ Rprintf("sampled phi\n"); }
    
    
    
    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postBranching") = postBranching,
                            Rcpp::Named("postBeta0") = postBeta0,
                            Rcpp::Named("postBeta") = postBeta,
                            Rcpp::Named("postTilde") = postTilde,
                            Rcpp::Named("postTau") = postTau,
                            Rcpp::Named("postZeta") = postZeta,
                            Rcpp::Named("postEta") = postEta,
                            Rcpp::Named("postPhi") = postPhi,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COVbeta0") = COVbeta0,
                            Rcpp::Named("COVbeta") = COVbeta,
                            Rcpp::Named("COVeta") = COVeta,
                            Rcpp::Named("COVphi") = COVphi);
}







// ESS for W, fixed rho, variance = 1, mean = 0
// MVN prior for beta0, beta, and delta
// Gamma prior for zeta and delta
// Unif prior for eta and phi
// [[Rcpp::export]]
List fitLGCPSE(int niter, vec ts, vec marks, List Xm, double maxT,
               double rho_beta, double rho_w,
               mat distmat, vec knts, mat tdiffm,
               double betaTilde0, vec beta0, vec betaTilde, mat beta, vec tau,
               vec Wm, double deltaTilde, vec delta, vec indlam0,
               vec zeta, double eta, double phi,
               double shape_tau, double rate_tau,
               double shape_zeta, double rate_zeta,
               double lb_eta, double ub_eta,
               double lb_phi, double ub_phi,
               vec sigma2, mat COVbeta0, List COVbeta,
               mat COVdelta, double COVeta, double COVphi,
               bool updateCOV, int adaptInterval, double adaptFactorExponent,
               vec adapIter, int ncores) {
  
  double negativeInf = -std::numeric_limits<float>::infinity();;
  // double positiveInf = std::numeric_limits<float>::infinity();;
  
  int n = ts.size(), p = beta.n_rows, K = distmat.n_rows, m = knts.size() - 1;
  double logprob;
  mat postBeta0 = zeros(niter, K), postWm(niter, (m+1)), postDelta(niter, K);
  mat postZeta = zeros(niter, K), postTau = zeros(niter, 1+p+1), postTilde = zeros(niter,p+2);
  vec postEta = zeros(niter), postPhi = zeros(niter);
  mat postBranching = zeros(niter, n), accprob = zeros(niter, 1+p+1+1+1);
  List postBeta(p);
  for(int i = 0; i < p; i ++){
    postBeta[i] = zeros(niter, K);
  }
  mat postBetai = zeros(niter, K);
  double sTilde, mTilde;
  vec postj;
  vec branching(n), numtrig(K);
  int count, numtrigtotal, numtrigl;
  uvec indtrig, indparent, indparentl;
  double tdifftrig, sdifftrig, dummySum, dummynewSum;
  int numback;
  vec lam0z0, newlam0z0, newWm(m+1);
  vec rhat = zeros(1+p+1+1+1), gamma1 = zeros(1+p+1+1+1), gamma2 = zeros(1+p+1+1+1), dummyaccprob;
  double c0 = 1, c1 = adaptFactorExponent, ropt = 0.234, shape, scale;
  vec xvals;
  double llprev, llthred, llnew, thetamin, thetamax, theta;
  vec priorWm, priorbeta0;
  bool accept;
  double zetal, neweta, newphi;
  vec newbeta0(K), betai(K), newbetai(K);
  mat newbeta(p, K);
  uvec indback, dummy;
  vec newdelta(K);
  
  mat cholCOVbeta0(K, K), cholCOVdelta(K, K);
  List cholCOVbeta(p);
  double cholCOVeta, cholCOVphi;
  cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
  for(int i = 0; i < p; i ++){
    cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
  }
  cholCOVdelta = trans( chol( sigma2[1+p] * ( COVdelta + 0.0000000001 * diagmat(ones(K)) ) ) );
  cholCOVeta = sqrt( sigma2[1+p+1] * (COVeta + 0.0000000001) );
  cholCOVphi = sqrt( sigma2[1+p+1+1] * (COVphi + 0.0000000001) );
  
  mat Vmat = exp(-distmat / rho_beta);
  mat Sigma = exp(-tdiffm / rho_w);
  mat invSigma = inv(Sigma), invVmat = inv(Vmat), cholSigma = chol(Sigma), cholVmat = chol(Vmat);
  double oneinvVmatone = sum( trans(ones(K)) * invVmat * ones(K) );
  vec invVmatone = invVmat * ones(K);
  
  double intLam0, newintLam0;
  vec lam0, newlam0;
  if(ncores == 1){
    intLam0 = compIntLam0(maxT, beta0, Xm, beta, exp(delta), Wm);
    lam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, exp(delta), Wm, indlam0, knts);
    
  } else {
    intLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, exp(delta), Wm, ncores);
    lam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, exp(delta), Wm, indlam0, knts, ncores);
  }
  
  
  
  
  // Start MCMC
  for(int s = 0; s < niter; s++) {
    
    // Gibbs update for latent branching structure
    branching = sampleBranching(ts, marks, distmat, lam0, zeta, eta, phi);
    postBranching.row(s) = trans( branching );
    if(s == 0){ Rprintf("sampled braching structure\n"); }
    
    
    // Self-exciting process
    indtrig = find( branching > 0 ); // triggered calls
    numtrigtotal = indtrig.size();
    indparent = indtrig;
    
    tdifftrig = 0;
    sdifftrig = 0;
    
    if(numtrigtotal != 0){
      for(int i = 0; i < numtrigtotal; i ++){
        indparent[i] = branching[ indtrig[i] ] - 1; // the delayed calls' parents
      }
      
      for(int l = 0; l < K; l ++){
        indparentl = find( (marks.elem(indparent) == l) );
        numtrigl = indparentl.size();
        
        if(numtrigl == 0){
          numtrig[l] = 0;
          
        } else {
          for(int i = 0; i < numtrigl; i ++){
            tdifftrig = tdifftrig + ts[ indtrig[ indparentl[i] ] ] - ts[ indparent[ indparentl[i] ] ];
            sdifftrig = sdifftrig + distmat(l, marks[ indtrig[ indparentl[i] ] ]);
          }
          numtrig[l] = numtrigl;
        }
      }
    }
    
    
    // Background process
    indback = find(branching == 0);
    lam0z0 = lam0.elem(indback);
    numback = indback.size();
    
    
    
    // M-H update for beta0
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(0);
        rhat[0] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[0] = 1 / pow(adapIter[0], c1);
        gamma2[0] = c0 * gamma1[0];
        sigma2[0] = exp( log(sigma2[0]) + gamma2[0] * (rhat[0] - ropt) );
        
        COVbeta0 = COVbeta0 + gamma1[0] * ( cov( postBeta0.rows(s+1-adaptInterval, s-1) ) - COVbeta0 );
        cholCOVbeta0 = trans( chol( sigma2[0] * ( COVbeta0 + 0.0000000001 * diagmat(ones(K)) ) ) );
        
        adapIter[0] = adapIter[0] + 1;
      }
    }
    
    newbeta0 = beta0 + cholCOVbeta0 * randn(K);
    if(ncores == 1){
      newlam0 = compLam0(ts, marks, maxT, newbeta0, Xm, beta, exp(delta), Wm, indlam0, knts);
      newintLam0 = compIntLam0(maxT, newbeta0, Xm, beta, exp(delta), Wm);
      
    } else {
      newlam0 = compLam0_parallel(ts, marks, maxT, newbeta0, Xm, beta, exp(delta), Wm, indlam0, knts, ncores);
      newintLam0 = compIntLam0_parallel(maxT, newbeta0, Xm, beta, exp(delta), Wm, ncores);
    }
    newlam0z0 = newlam0.elem(indback);
    
    logprob = (sum(log(newlam0z0)) - sum(log(lam0z0))) - (newintLam0 - intLam0) +
      MVN_logh(newbeta0, betaTilde0 * ones(K), invVmat/tau[0]) - MVN_logh(beta0, betaTilde0 * ones(K), invVmat/tau[0]);
    
    if (log(randu()) < logprob) {
      beta0 = newbeta0;
      lam0 = newlam0;
      intLam0 = newintLam0;
      lam0z0 = newlam0z0;
      
      accprob(s, 0) = 1;
    }
    
    postBeta0.row(s) = trans(beta0);
    if(s == 0){ Rprintf("sampled beta0\n"); }
    
    
    
    // Gibbs update for betaTilde0
    sTilde = 1 / (oneinvVmatone / tau[0] + 1/100);
    mTilde = sTilde / tau[0] * sum(trans(beta0) * invVmatone);
    betaTilde0 = mTilde + sqrt(sTilde) * randn();
    postTilde(s,0) = betaTilde0;
    if(s == 0){ Rprintf("sampled betaTilde0\n"); }
    
    
    
    // M-H updates for beta
    for(int i = 0; i < p; i ++){
      if( updateCOV ){
        if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
          dummyaccprob = accprob.col(1+i);
          rhat[1+i] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
          gamma1[1+i] = 1 / pow(adapIter[1+i], c1);
          gamma2[1+i] = c0 * gamma1[1+i];
          sigma2[1+i] = exp( log(sigma2[1+i]) + gamma2[1+i] * (rhat[1+i] - ropt) );
          
          postBetai = as<mat>(postBeta[i]);
          COVbeta[i] = as<mat>(COVbeta[i]) + gamma1[1+i] * ( cov( postBetai.rows(s+1-adaptInterval, s-1) ) - as<mat>(COVbeta[i]) );
          cholCOVbeta[i] = trans( chol( sigma2[1+i] * ( as<mat>(COVbeta[i]) + 0.0000000001 * diagmat(ones(K)) ) ) );
          
          adapIter[1+i] = adapIter[1+i] + 1;
        }
      }
      
      betai = trans(beta.row(i));
      newbetai = betai + as<mat>(cholCOVbeta[i]) * randn(K);
      newbeta = beta;
      newbeta.row(i) = trans(newbetai);
      
      if(ncores == 1){
        newlam0 = compLam0(ts, marks, maxT, beta0, Xm, newbeta, exp(delta), Wm, indlam0, knts);
        newintLam0 = compIntLam0(maxT, beta0, Xm, newbeta, exp(delta), Wm);
        
      } else {
        newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, newbeta, exp(delta), Wm, indlam0, knts, ncores);
        newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, newbeta, exp(delta), Wm, ncores);
      }
      newlam0z0 = newlam0.elem(indback);
      
      logprob = (sum(log(newlam0z0)) - sum(log(lam0z0))) - (newintLam0 - intLam0) +
        MVN_logh(newbetai, betaTilde[i] * zeros(K), invVmat / tau[1+i]) - MVN_logh(betai, betaTilde[i] * zeros(K), invVmat / tau[1+i]);
      
      if (log(randu()) < logprob) {
        betai = newbetai;
        beta = newbeta;
        lam0 = newlam0;
        intLam0 = newintLam0;
        lam0z0 = newlam0z0;
        
        accprob(s, 1+i) = 1;
      }
      
      postBetai = as<mat>(postBeta[i]);
      postBetai.row(s) = trans(betai);
      postBeta[i] = postBetai;
    }
    if(s == 0){ Rprintf("sampled beta\n"); }
    
    
    
    // Gibbs update for betaTilde
    for(int i = 0; i < p; i ++){
      sTilde = 1 / (oneinvVmatone / tau[1+i] + 1/100);
      mTilde = sTilde / tau[1+i] * sum(beta.row(i) * invVmatone);
      betaTilde[i] = mTilde + sqrt(sTilde) * randn();
      postTilde(s,1+i) = betaTilde[i]; 
    }
    if(s == 0){ Rprintf("sampled betaTilde\n"); }
    
    
    
    // Elliptical slice sampling for Wm
    llprev = sum(log(lam0z0)) - intLam0;
    priorWm = trans( cholSigma ) * randn(m+1);
    
    thetamin = 0;
    thetamax = 2 * M_PI;
    theta = thetamin + randu() * (thetamax - thetamin);
    thetamin = theta - 2 * M_PI;
    thetamax = theta;
    
    llthred = llprev + log(randu());
    accept = false;
    count = 0;
    
    while(accept == false){
      count = count + 1;
      
      newWm = Wm * cos(theta) + priorWm * sin(theta);
      
      if(ncores == 1){
        newlam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, exp(delta), newWm, indlam0, knts);
        newintLam0 = compIntLam0(maxT, beta0, Xm, beta, exp(delta), newWm);
        
      } else {
        newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, exp(delta), newWm, indlam0, knts, ncores);
        newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, exp(delta), newWm, ncores);
      }
      
      newlam0z0 = newlam0.elem(indback);
      llnew = sum(log(newlam0z0)) - newintLam0;
      
      if(llnew > llthred){
        llprev = llnew;
        accept = true;
      } else {
        if(theta < 0){ thetamin = theta; } else { thetamax = theta; }
        theta = thetamin + randu() * (thetamax - thetamin);
        if ( (count) % 1000 == 0 ) {
          Rprintf("ESS for Wm: %d iterations...\n", count);
        }
      }
    }
    Wm = newWm;
    lam0 = newlam0;
    intLam0 = newintLam0;
    lam0z0 = newlam0z0;
    
    postWm.row(s) = trans(Wm);
    if(s == 0){ Rprintf("sampled W\n"); }
    
    
    
    // MH updates for delta
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1+p);
        rhat[1+p] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1+p] = 1 / pow(adapIter[1+p], c1);
        gamma2[1+p] = c0 * gamma1[1+p];
        sigma2[1+p] = exp( log(sigma2[1+p]) + gamma2[1+p] * (rhat[1+p] - ropt) );
        
        COVdelta = COVdelta + gamma1[1+p] * ( cov( postDelta.rows(s+1-adaptInterval, s-1) ) - COVdelta );
        cholCOVdelta = trans( chol( sigma2[1+p] * ( COVdelta + 0.0000000001 * diagmat(ones(K)) ) ) );
        
        adapIter[1+p] = adapIter[1+p] + 1;
      }
    }
    
    newdelta = delta + cholCOVdelta * randn(K);
    
    if(ncores == 1){
      newlam0 = compLam0(ts, marks, maxT, beta0, Xm, beta, exp(newdelta), Wm, indlam0, knts);
      newintLam0 = compIntLam0(maxT, beta0, Xm, beta, exp(newdelta), Wm);
      
    } else {
      newlam0 = compLam0_parallel(ts, marks, maxT, beta0, Xm, beta, exp(newdelta), Wm, indlam0, knts, ncores);
      newintLam0 = compIntLam0_parallel(maxT, beta0, Xm, beta, exp(newdelta), Wm, ncores);
    }
    newlam0z0 = newlam0.elem(indback);
    
    logprob = (sum(log(newlam0z0)) - sum(log(lam0z0))) - (newintLam0 - intLam0) +
      MVN_logh(newdelta, deltaTilde * ones(K), invVmat/tau[1+p]) - MVN_logh(delta, deltaTilde * ones(K), invVmat/tau[1+p]);
    
    if (log(randu()) < logprob) {
      delta = newdelta;
      lam0 = newlam0;
      intLam0 = newintLam0;
      lam0z0 = newlam0z0;
      
      accprob(s, 1+p) = 1;
    }
    
    postDelta.row(s) = trans(delta);
    if(s == 0){ Rprintf("sampled delta\n"); }
    
    
    
    // Gibbs update for deltaTilde
    sTilde = 1 / (oneinvVmatone / tau[1+p] + 1/100);
    mTilde = sTilde / tau[1+p] * sum(trans(delta) * invVmatone);
    deltaTilde = mTilde + sqrt(sTilde) * randn();
    postTilde(s,1+p) = deltaTilde;
    if(s == 0){ Rprintf("sampled deltaTilde\n"); }
    
    
    
    // Gibbs update for tau
    shape = shape_tau + 0.5 * K;
    scale = ( 1 / ( rate_tau + 0.5 * trans(beta0 - betaTilde0) * invVmat * (beta0 - betaTilde0) ) )[0];
    tau[0] = 1/randg(distr_param(shape, scale));
    
    for(int i = 0; i < p; i ++){
      shape = shape_tau + 0.5 * K;
      scale = ( 1 / ( rate_tau + 0.5 * (beta.row(i) - betaTilde[i]) * invVmat * trans(beta.row(i) - betaTilde[i]) ) )[0];
      tau[1+i] = 1/randg(distr_param(shape, scale));
    }
    
    shape = shape_tau + 0.5 * K;
    scale = ( 1 / ( rate_tau + 0.5 * trans(delta - deltaTilde) * invVmat * (delta - deltaTilde) ) )[0];
    tau[1+p] = 1/randg(distr_param(shape, scale));
    
    postTau.row(s) = trans(tau);
    if(s == 0){ Rprintf("sampled tau\n"); }
    
    
    
    // Gibbs update for zeta
    for(int l = 0; l < K; l ++){
      
      if(ncores == 1){
        dummySum = compSumIntHl(l, ts, marks, maxT, distmat, eta, phi);
      } else{
        dummySum = compSumIntHl_parallel(l, ts, marks, maxT, distmat, eta, phi, ncores);
      }
      
      shape = shape_zeta + numtrig[l];
      scale = 1 / ( rate_zeta +  dummySum );
      
      zetal = randg(distr_param(shape, scale));
      zeta[l] = zetal;
      
      postZeta(s,l) = zetal;
    }
    if(s == 0){ Rprintf("sampled zeta\n"); }
    
    
    
    // M-H update for eta
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1+p+1);
        rhat[1+p+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1+p+1] = 1 / pow(adapIter[1+p+1], c1);
        gamma2[1+p+1] = c0 * gamma1[1+p+1];
        sigma2[1+p+1] = exp( log(sigma2[1+p+1]) + gamma2[1+p+1] * (rhat[1+p+1] - ropt) );
        
        COVeta = COVeta + gamma1[1+p+1] * ( var( postEta.rows(s+1-adaptInterval, s-1) ) - COVeta );
        cholCOVeta = sqrt( sigma2[1+p+1] * (COVeta + 0.0000000001) );
        
        adapIter[1+p+1] = adapIter[1+p+1] + 1;
      }
    }
    
    neweta = eta + cholCOVeta * randn();
    
    if( (neweta < lb_eta) || (neweta > ub_eta) ){
      logprob = negativeInf;
      
    } else {
      if(ncores == 1){
        dummynewSum = compSumIntH(ts, marks, maxT, distmat, zeta, neweta, phi);
        dummySum = compSumIntH(ts, marks, maxT, distmat, zeta, eta, phi);
      } else{
        dummynewSum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, neweta, phi, ncores);
        dummySum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, eta, phi, ncores);
      }
      
      if( numtrigtotal == 0 ){
        logprob = - ( dummynewSum - dummySum );
        
      } else {
        logprob = - tdifftrig * (neweta - eta) - ( dummynewSum - dummySum );
      }
    }
    
    if (log(randu()) < logprob) {
      eta = neweta;
      accprob(s, 1+p+1) = 1;
    }
    
    postEta[s] = eta;
    if(s == 0){ Rprintf("sampled eta\n"); }
    
    
    
    // M-H update for phi
    if( updateCOV ){
      if( (s+1 >= adaptInterval) && (s+1 - (adaptInterval * trunc((s+1) / adaptInterval)) == 0) ){
        dummyaccprob = accprob.col(1+p+1+1);
        rhat[1+p+1+1] = sum( dummyaccprob.rows(s+1-adaptInterval, s-1) ) / (adaptInterval-1);
        gamma1[1+p+1+1] = 1 / pow(adapIter[1+p+1+1], c1);
        gamma2[1+p+1+1] = c0 * gamma1[1+p+1+1];
        sigma2[1+p+1+1] = exp( log(sigma2[1+p+1+1]) + gamma2[1+p+1+1] * (rhat[1+p+1+1] - ropt) );
        
        COVphi = COVphi + gamma1[1+p+1+1] * ( var( postPhi.rows(s+1-adaptInterval, s-1) ) - COVphi );
        cholCOVphi = sqrt( sigma2[1+p+1+1] * (COVphi + 0.0000000001) );
        
        adapIter[1+p+1+1] = adapIter[1+p+1+1] + 1;
      }
    }
    
    newphi = phi + cholCOVphi * randn();
    
    if( (newphi < lb_phi) || (newphi > ub_phi) ){
      logprob = negativeInf;
      
    } else {
      if(ncores == 1){
        dummynewSum = compSumIntH(ts, marks, maxT, distmat, zeta, eta, newphi);
        dummySum = compSumIntH(ts, marks, maxT, distmat, zeta, eta, phi);
      } else{
        dummynewSum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, eta, newphi, ncores);
        dummySum = compSumIntH_parallel(ts, marks, maxT, distmat, zeta, eta, phi, ncores);
      }
      
      if( numtrigtotal == 0 ){
        logprob = - ( dummynewSum - dummySum );
        
      } else {
        logprob = - sdifftrig * (newphi - phi) - ( dummynewSum - dummySum  );
      }
    }
    
    if (log(randu()) < logprob) {
      phi = newphi;
      accprob(s, 1+p+1+1) = 1;
    }
    
    postPhi[s] = phi;
    if(s == 0){ Rprintf("sampled phi\n"); }
    
    
    
    if ( (s+1) % 10 == 0 ) {
      Rprintf("Generated %d samples...\n", s+1);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("postBranching") = postBranching,
                            Rcpp::Named("postBeta0") = postBeta0,
                            Rcpp::Named("postBeta") = postBeta,
                            Rcpp::Named("postTau") = postTau,
                            Rcpp::Named("postDelta") = postDelta,
                            Rcpp::Named("postTilde") = postTilde,
                            Rcpp::Named("postWm") = postWm,
                            Rcpp::Named("postZeta") = postZeta,
                            Rcpp::Named("postEta") = postEta,
                            Rcpp::Named("postPhi") = postPhi,
                            Rcpp::Named("Wm") = Wm,
                            Rcpp::Named("Accprob") = accprob,
                            Rcpp::Named("sigma2") = sigma2,
                            Rcpp::Named("adapIter") = adapIter,
                            Rcpp::Named("COVbeta0") = COVbeta0,
                            Rcpp::Named("COVbeta") = COVbeta,
                            Rcpp::Named("COVdelta") = COVdelta,
                            Rcpp::Named("COVeta") = COVeta,
                            Rcpp::Named("COVphi") = COVphi);
}



