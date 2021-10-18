mu_d_beta <- function(i,
                      j,
                      q,
                      gammas,
                      B,
                      X,
                      Z,
                      P,
                      X_tilde,
                      Z_tilde,
                      Z_tilde_gamma_cols,
                      P_tilde,
                      gamma_tilde){


  mu_deriv <-
    t(X[i,q,drop = F])%*%(Z[i,,drop = F]%*%P[,j,drop = F]*exp(gammas[i] +
                                                                X[i,,drop = F]%*%B[,j,drop = F]))

  K_tilde <- dim(P_tilde)[1]
  for(k_tilde in 1:K_tilde){
    if(k_tilde %in% Z_tilde_gamma_cols){
      mu_deriv <- mu_deriv + exp(gammas[i])*
        (Z_tilde[i,k_tilde,drop = F])%*%
        P_tilde[k_tilde,j,drop = F]*
        exp(gamma_tilde[k_tilde] +
              X_tilde[k_tilde,,drop = F]%*%B[,j,drop = F])*
        X_tilde[k_tilde,q,drop = F]
    } else{
      mu_deriv <- mu_deriv +
        (Z_tilde[i,k_tilde,drop = F])%*%
        P_tilde[k_tilde,j,drop = F]*
        exp(gamma_tilde[k_tilde] +
              X_tilde[k_tilde,,drop = F]%*%B[,j,drop = F])*
        X_tilde[k_tilde,q,drop = F]
    }
  }



  return(mu_deriv)
}
