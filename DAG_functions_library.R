#### packages needed ####
library(CholWishart)
library(BCDAG)
library(mvtnorm)
library(graph)
library(BiocGenerics)
library(gRbase)
library(Rgraphviz)

### sources ####
source("propose_DAG.R")
source("operation.R")
source("pa.R")
source("fa.R")

##Functions #####

Gamma_f = function(p, a){
  gama <- rep(NA,p)
  if(p == 0){
    return(0)
  } else{
    for(i in 1:p){
      gama[i] <- lgamma(a+(1-i)*0.5)
    }
    gamma_out <- 0.25*p*(p-1)*log(pi) + sum(gama)
    return(gamma_out)
  }
}


PosteriorDraws <- function(Y,X, B0, C, U){
  
  a = q
  
  n = length(Y[,1])
  
  B_hat <- solve(t(X)%*%X)%*%t(X)%*%Y
  
  B_tilde <- solve(t(X)%*%X + C)%*%(t(X)%*%Y + C%*%B0)
  
  E_hat <- Y - X%*%B_hat
  
  D <- t(B_hat)%*%t(X)%*%X%*%B_hat + t(B0)%*%C%*%B0 - t(C%*%B0 + t(X)%*%Y)%*%solve(t(X)%*%X + C)%*%(C%*%B0 + t(X)%*%Y)
  
  
  C_tilde <- t(X)%*%X + C
  
  
  U_tilde = U + t(E_hat)%*%E_hat + D
  
  a_tilde = a + n
  
  
  Omega_post = rWishart(1, a_tilde, solve(U_tilde))[,,1]
  
  B_tilde_vec = c(B_tilde)
  
  Omega_C_tilde = kronecker(solve(Omega_post), solve(C_tilde)) 
  
  B_post = matrix(rmvnorm(1, B_tilde_vec, Omega_C_tilde), p, q)
  
  
  return(out_post = list(B_post = B_post, Omega_post = Omega_post))
  
}


Norm_const <- function(C, U, a, q){
  q*(p+1)/2*log(2*pi)+(a*q/2)*log(2)+Gamma_f(q,a/2)-
    (q/2)*log(det(C))-(a/2)*log(det(U))
}


marginal_lik = function(A, Y, Z, a, U, B0, C){
  if(length(A)== 0){marg = 0}else{
    
    Y_A = Y[,A,drop = FALSE]
    U_A = U[A,A,drop = FALSE]
    
    A_a = (q - length(A))
    n = nrow(Y) #domanda per prof
    
    B0_A = B0[,A,drop = FALSE]
    
    B_hat <- solve(t(Z)%*%Z)%*%t(Z)%*%Y_A
    B_tilde <- solve(t(Z)%*%Z + C)%*%(t(Z)%*%Y_A + C%*%B0_A)
    E_hat <- Y_A - Z%*%B_hat
    D <- t(B_hat)%*%t(Z)%*%Z%*%B_hat + t(B0_A)%*%C%*%B0_A - t(C%*%B0_A + t(Z)%*%Y_A)%*%solve(t(Z)%*%Z + C)%*%(C%*%B0_A + t(Z)%*%Y_A)
    C_tilde <- t(Z)%*%Z + C
    U_tilde = U_A + t(E_hat)%*%E_hat + D
    A_tilde = A_a + n
    p_tilde = nrow(Y) # check
    
    priorc = Norm_const(C = C, U = U_A, a = a - (q - length(A)), q = length(A))
    
    posteriorc = Norm_const(C_tilde, U_tilde, A_tilde, p_tilde)
    
    marg = posteriorc - priorc - (n*q)/2*log(2*pi)
    
  }
  
  return(marg)
  
}


marginal_lik_revised = function(A, a, q, n, U, C, U_tilde, C_tilde){
  
  if(length(A)== 0){marg = 0}else{
    
    U_A = U[A,A,drop = FALSE]
    U_tilde_A = U_tilde[A,A,drop = FALSE]
    
    a_A = a - (q - length(A))
    
    priorc = Norm_const(C = C, U = U_A, a = a_A, q = length(A))
    
    posteriorc = Norm_const(C = C_tilde, U = U_tilde_A, a = a_A + n, q = length(A))
    
    marg = posteriorc - priorc - (n*length(A))/2*log(2*pi)
    
  }
  
  return(marg)
  
}


marg_DAG <- function(DAG,Y,Z){
  sum(sapply(1:q, function(j) marginal_lik(A = fa(j, DAG), Y = Y, Z = Z, a = a, U = U, B0 = B0, C = C))) -
    sum(sapply(1:q, function(j) marginal_lik(A = pa(j, DAG), Y = Y, Z = Z, a = a, U = U, B0 = B0, C = C)))
}

marg_DAG_revised <- function(DAG, a, q, n, U, C, U_tilde, C_tilde){
  sum(sapply(1:q, function(j) marginal_lik_revised(A = fa(j, DAG), a, q, n, U, C, U_tilde, C_tilde))) -
    sum(sapply(1:q, function(j) marginal_lik_revised(A = pa(j, DAG), a, q, n, U, C, U_tilde, C_tilde)))
}

mcmc_dag_revised = function(Y, Z, S, w){
  
  # Y, Z : data
  # S    : number of MCMC iterations
  # w    : prior probability of edge inclusion in p(D)
  
  # plus other parameters if needed
  
  # Require functions to compute marginal likelihood of a DAG, perform moves between DAGs (proposal), ...
  # External function written in R file "script.R" can be loaded using command source("script.R")
  
  # For instance:
  
  # souce(marg_like_dag.R) (once implemented!)
  
  n = dim(Y)[1]
  q = dim(Y)[2]
  
  # a = q + 1
  
  a = q
  
  Graph_post = array(NA, c(q, q, S)) # here the S adjacency matrices of visited DAGs will be stored
  
  num_edges = c()
  
  # Initialize the chain at the empty DAG
  
  DAG = matrix(0, q, q)
  
  Graph_post[,,1] = DAG
  
  ########################
  ## Set hyperparamters ##
  ########################
  
  q = ncol(Y)
  p = ncol(Z)
  kq = q
  
  a = q
  U = diag(0.01, q)
  B0 = matrix(0, p, q)
  C = diag(0.01, p)
  
  ###################################
  ## Compute sufficient statistics ##
  ###################################
  
  B_hat <- solve(t(Z)%*%Z)%*%t(Z)%*%Y
  B_tilde <- solve(t(Z)%*%Z + C)%*%(t(Z)%*%Y + C%*%B0)
  E_hat <- Y - Z%*%B_hat
  
  #D <- t(B_hat)%*%t(Z)%*%Z%*%B_hat + t(B0)%*%C%*%B0 - t(C%*%B0 + t(Z)%*%Y)%*%solve(t(Z)%*%Z + C)%*%(C%*%B0 + t(Z)%*%Y)
  
  D <- t(B0 - B_hat)%*%solve(solve(C) + solve(t(Z)%*%Z))%*%(B0 - B_hat)
  
  C_tilde <- t(Z)%*%Z + C
  U_tilde = U + t(E_hat)%*%E_hat + D
  
  for(s in 1:S){
    
    ## Update the graph
    
    Graph_move = propose_DAG(DAG, fast = TRUE) # propose a new DAG using the function
    
    DAG_star = Graph_move$proposedDAG
    
    # Evaluate the priors p(D_star) and p(DAG) (equivalently the logs) if not uniform (otherwise set 0)
    
    prior_DAG = sum(DAG)*log(w) + (q*(q-1)/2 - sum(DAG))*log(1-w)
    prior_DAG_star = sum(DAG_star)*log(w) + (q*(q-1)/2 - sum(DAG_star))*log(1-w)
  
    # prior_DAG_star = 0
    # prior_DAG      = 0
    
    # Includere prior su DAG
    
    # Compute the marginal likelihood of DAG_star and DAG with the function marg_like_dag (una volta costruita la funzione!)
    
    margin_star = marg_DAG_revised(DAG = DAG_star, a, q, n, U, C, U_tilde, C_tilde)
    margin      = marg_DAG_revised(DAG = DAG, a, q, n, U, C, U_tilde, C_tilde)
    
    # Compute log acceptance ratio
    
    ratio_D = min(0, margin_star - margin + prior_DAG_star - prior_DAG)
    
    # accept the proposed DAG
    
    if(log(runif(1)) < ratio_D){
      
      DAG = DAG_star
      
    }
    
    # Store the DAG in the output
    
    Graph_post[,,s] = DAG
    
    num_edges[s] = sum(DAG)
    
  }
  
  return(out = list(Y = Y,
                    Z = Z,
                    Graph_post = Graph_post,
                    num_edges = num_edges))
  
}
