data {
  int<lower=1> n;//
  int<lower=1> K;//
  matrix[K,n] X;//
  matrix[n,n] G;//
  vector<lower=0>[n] n_p;//
  int<lower=1> topo[n];//
  int<lower=0> topo_h[n];//
  int<lower=0> topo_2[n];//
  int<lower=0> topo_1[n];//
  int<lower=0> topo_no1[n];//
  real cauchy_par;//
}

parameters {
  real<lower=0> tau;//the logarithm of the global scale parameter
  real<lower=0, upper=1> psi;//decay factor
  //matrix[K,n] delta;//delta_ij
  vector[n] eta;//
  matrix[K,n] w;//
  // p=1 
  real<lower=0> sigma2;
}

transformed parameters{
  vector[n] phi;
  vector[n] h;//
  matrix[K,n] Gdiff_1;//
  matrix[K,n] Gdiff_2;//
  matrix[K,n] mu;//
  
  //h
  //
   for(i in 1:n){
    h[i]=2*log(tau)+eta[i];
    phi[i]=0;
  }
  
  //
  for(i in topo_h){
    if(i!=0){
      phi[i]=psi*(dot_product(G[,i],h)/n_p[i]-2*log(tau));
      h[i]+=phi[i];
    }
  }
  
  
  //
  Gdiff_2=rep_matrix(0,K,n);
  for(i in topo_h){
    if(i!= 0){
      Gdiff_2[,i]=w[,i];
    }
  }
 
  //
  Gdiff_1=Gdiff_2;
  for(i in topo_no1){
    if(i!= 0){
      Gdiff_1[,i]=w[,i];
    }
  }
   for(i in topo_h){
    if(i!=0){
        Gdiff_1[,i]=Gdiff_2[,i]+(Gdiff_1*G[,i])/n_p[i];
    }
  }
  //mu
  //
  mu=rep_matrix(0,K,n);
  for(i in topo_1){
    if(i!=0){
        mu[,i]=w[,i];
    }
  }
  //
  for(i in topo_no1){
    if(i!=0){
      mu[,i]=Gdiff_1[,i]+(mu*G[,i])/n_p[i];
    }
  }

}

model {
  //
  target+=eta/2-log1p_exp(eta);
  //
  psi~uniform(0,1);
  //tau
  tau~cauchy(0,cauchy_par);
  
  //delta

  
  //w~normal
  for(i in 1:n){
    for(j in 1:K){
       target+=normal_lpdf(w[j,i] | 0,sqrt(exp(h[i])));
    }
  }
  
  //sigma
  // p=1
  sigma2 ~ inv_gamma(1e-2, 1e-2);
  
  //
  // p=1
  for(i in 1:n){
    for(j in 1:K){
      target += normal_lpdf(X[j,i] | mu[j,i],sqrt(sigma2));
    }
  }
 
}
