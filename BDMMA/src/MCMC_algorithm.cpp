#include <iostream>
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <math.h>



//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppEigen)]]


using namespace Rcpp;


//[[Rcpp::export]]

NumericVector DM(arma::mat Y, arma::mat alpha){
  NumericVector logl;
  arma::vec alpha_rowsums,m;
  alpha_rowsums=sum(alpha,1);
  m=sum(Y,1);
  logl = sum(lgamma(alpha_rowsums) - lgamma(alpha_rowsums+m))
    + sum(sum(lgamma(alpha+Y),1)) - sum(sum(lgamma(alpha),1));

  return(logl);
}



//[[Rcpp::export]]

double DM_diff(arma::mat Y, arma::mat alpha1, arma::mat alpha2, int g){
  double logl;
  arma::vec alpha_rowsums1,alpha_rowsums2,m;
  alpha_rowsums1=sum(alpha1,1);
  alpha_rowsums2=sum(alpha2,1);
  m=sum(Y,1);
  logl = sum(lgamma(alpha_rowsums1) - lgamma(alpha_rowsums2) - lgamma(alpha_rowsums1+m) +
    lgamma(alpha_rowsums2 + m) + lgamma(alpha1.col(g)+Y.col(g)) - lgamma(alpha2.col(g)+Y.col(g)) -
    lgamma(alpha1.col(g)) + lgamma(alpha2.col(g)));
  return(logl);
}



//[[Rcpp::export]]

arma::mat update_alpha2(arma::rowvec alpha, arma::mat alpha_m,arma::mat x,
                        arma::mat y, arma::mat beta, int T, double lambda,
                        arma::vec prop, arma::mat delta_m){
  NumericVector f;
  double alpha_g;
  arma::mat alpha_m2 = alpha_m;
  arma::mat gamma1, gamma2;
  for (int g=0;g<T;g++){
    alpha_g=exp(log(alpha(g))+R::rnorm(0,0.02));
    alpha_m2=alpha_m;
    alpha_m2.col(g)=alpha_m2.col(g)+(alpha_g-alpha(g));
    gamma1=alpha_m % exp(x * beta + delta_m);
    gamma2=alpha_m2 % exp(x * beta + delta_m);
    f=exp(DM_diff(y,gamma2,gamma1,g)+(alpha[g]-alpha_g)/(lambda*prop[g]));
    if (is_true(any(runif(1,0,1)<f))){
      alpha[g]=alpha_g;
      alpha_m=alpha_m2;
    }
  }
  return (alpha_m);
}



//[[Rcpp::export]]

arma::mat update_beta1(arma::mat alpha_m, arma::mat x, arma::mat y,
                       arma::mat beta, int T, arma::vec L, double sigma1,
                       double eta, arma::mat delta_m){
  double f;
  double beta_g, f1;
  arma::mat beta2 = beta;
  arma::mat gamma1, gamma2;
  for (int g=0;g<T;g++){
    beta_g=beta(0,g)+R::rnorm(0,0.02);
    beta2(0,g)=beta_g;
    gamma1=alpha_m % exp(x * beta + delta_m);
    gamma2=alpha_m % exp(x * beta2 + delta_m);
    f1=DM_diff(y,gamma2,gamma1,g);
    f=f1+L(g)*(-pow(beta2(0,g),2)/(2*pow(sigma1,2))+pow(beta(0,g),2)/(2*pow(sigma1,2)))+
      (1-L(g))*(-pow(beta2(0,g),2)+pow(beta(0,g),2))/(pow(eta,2)*2);
    if (is_true(any(runif(1,0,1)<exp(f)))){
      beta(0,g)=beta_g;
    }
    beta2=beta;
  }
  return(beta);
}



//[[Rcpp::export]]

arma::mat update_beta(arma::mat alpha_m, arma::mat x, arma::mat y,
                      arma::mat beta, int T, int K, double sigma2, arma::mat delta_m){
  double f;
  double beta_g,f1;
  arma::mat beta2 =  beta;
  arma::mat gamma1, gamma2;
  for (int g=0;g<T;g++){
    for (int j=1;j<K;j++){
      beta_g=beta(j,g)+R::rnorm(0,0.02);
      beta2(j,g)=beta_g;
      gamma1=alpha_m % exp(x * beta + delta_m);
      gamma2=alpha_m % exp(x * beta2 + delta_m);
      f1=DM_diff(y,gamma2,gamma1,g);
      f=f1+(pow(beta(j,g),2)-pow(beta2(j,g),2))/(2*pow(sigma2,2));
      if (is_true(any(runif(1,0,1)<exp(f)))){
        beta(j,g)=beta_g;
      }
      beta2=beta;
    }
  }
  return(beta);
}




//[[Rcpp::export]]

arma::mat update_delta(arma::mat alpha_m, arma::mat x, arma::mat y, arma::mat delta,
                       arma::mat delta_m, arma::mat beta, arma::mat e_delta,
                       arma::vec batch, int N, int T, int I, double sigma3){
  NumericVector f;
  double delta_g;
  arma::mat delta2 = delta;
  arma::mat gamma1, gamma2, delta_m2 = delta_m;
  for (int g=0;g<T;g++){
    for (int j=0;j<(I-1);j++){
      delta_g=delta(j,g)+R::rnorm(0,0.02);
      delta2(j,g)=delta_g;
      for (int i=0;i<(N-1);i++) {
        delta_m2.row(i)=delta2.row(batch(i+1)-1);
      }
      gamma1=alpha_m % exp(x * beta + delta_m);
      gamma2=alpha_m % exp(x * beta + delta_m2);
      f=DM_diff(y,gamma2,gamma1,g)-(pow((delta2(j,g)-e_delta(g,j)),2)-pow((delta(j,g)-e_delta(g,j)),2))/(2*pow(sigma3,2));
      if (is_true(any(runif(1,0,1)<exp(f)))){
        delta(j,g)=delta_g;
      }
      delta.row(3)=-(109*delta.row(0)+152*delta.row(1)+165*delta.row(2))/100;
      for (int i=0;i<(N-1);i++) {
        delta_m.row(i)=delta.row(batch(i+1)-1);
      }
      delta_m2=delta_m;
      delta2=delta;
    }
  }
  return(join_cols(delta,delta_m));
}



//[[Rcpp::export]]

arma::rowvec Mat_To_Rowvec(arma::mat m){
  arma::vec v(m.n_cols*m.n_rows);
  v=arma::vectorise(m);
  return(v.t());
}


//[[Rcpp::export]]

arma::mat MCMC(arma::rowvec alpha, arma::mat alpha_m, arma::mat x, arma::mat y,
              arma::mat beta, arma::mat delta, arma::mat delta_m, arma::mat e_delta,
              int T, int N, int K, int I, double lambda, arma::vec prop, arma::vec L,
              double sigma1, double sigma2, double sigma3, int iter,
              double eta, double a, double b, double p, arma::vec batch){

  double shape, scale, shape2, scale2, prob, p1, p2, var;
  arma::mat delta_out, final_out, alpha_all(iter,T), beta1_all(iter,T), beta2_all(iter,T*(K-1)),
  L_all(iter,T), beta_out(K+1,T), delta_all(iter,T*I);
  arma::rowvec mean_alpha(T), mean_beta1(T);
  arma::vec mean_L(T), mean_p(T), eta_all(iter), p_all(iter);

  for (int iteration = 1; iteration <= iter; iteration++){

    // Update alpha
    alpha_m = update_alpha2(alpha = alpha,alpha_m = alpha_m,x = x,y = y,
                            beta=beta, T=T, lambda = lambda, prop = prop,
                            delta_m=delta_m);
    alpha=alpha_m.row(0);
    alpha_all.row(iteration-1)=alpha;

    // Update Lambda
    shape=T-1;
    scale=(1/sum(alpha.t()/prop));
    lambda=1/(R::rgamma(shape,scale));

    // Update beta[1,]
    beta=update_beta1(alpha_m = alpha_m,x=x,y=y,beta = beta,T=T,L=L,
                      sigma1=sigma1,eta=eta,delta_m=delta_m);
    beta1_all.row(iteration-1)=beta.row(0);

    // Update L
    for (int g=0;g<T;g++){
      p1=p/sigma1*exp(-pow(beta(0,g),2)/(2*pow(sigma1,2)));
      p2=(1-p)/eta*exp(-pow(beta(0,g),2)/(2*pow(eta,2)));
      prob=p1/(p1+p2);
      L(g)=R::rbinom(1,prob);
    }

    L_all.row(iteration-1)=L.t();

    //Update p
    p=R::rbeta(a+sum(L),b+sum(1-L));
    p_all(iteration-1)=p;


    // Update eta
    shape2=0.1+sum(1-L)/2;
    scale2=0.1+sum((1-L) % pow((beta.row(0)).t(),2))/2;
    eta=pow(1/(R::rgamma(shape2,scale2)),0.5);
    eta_all(iteration-1)=eta;

    // Update beta[-1,]
    if (K>=2){
      beta=update_beta(alpha_m = alpha_m,x=x,y=y,beta = beta,T=T,
                       K=K,sigma2 = sigma2,delta_m=delta_m);
      beta2_all.row(iteration-1)=Mat_To_Rowvec(beta.rows(1,K-1));
    }

    // Update delta
    delta_out=update_delta(alpha_m=alpha_m, x=x, y=y, delta=delta,
                           delta_m=delta_m, beta=beta, e_delta=e_delta,
                           batch=batch, N=N, T=T, I=I, sigma3=sigma3);

    delta=delta_out.rows(0,I-1);
    delta_m=delta_out.rows(I,N+I-1);
    delta_all.row(iteration-1)=arma::vectorise(delta).t();

    if (iteration % 50 == 0){
      Rcout << "Iteration = " << iteration <<"\n";
    }
  }
  final_out=join_rows(alpha_all,beta1_all);
  final_out=join_rows(final_out,L_all);
  final_out=join_rows(join_rows(final_out,eta_all),p_all);
  final_out=join_rows(final_out,beta2_all);
  return(join_rows(final_out,delta_all));
}

