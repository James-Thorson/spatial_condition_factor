// Space time 
#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Dimensions
  DATA_INTEGER(n_i);         // Total number of observations
  DATA_INTEGER(n_s);         // Number of "strata" (i.e., vectices in SPDE mesh) 
  DATA_INTEGER(n_x)	        // Number of stations (i.e., in k-means algorithm) 
  DATA_INTEGER(n_t)          // Number of years  
  DATA_INTEGER(n_j_a)	        // Number of covariates
  DATA_INTEGER(n_j_b)	        // Number of covariates

  // Config
  DATA_INTEGER(Aniso); // 0: No; 1:Yes
  DATA_FACTOR(ObsModel);    // Observation model
  DATA_FACTOR(FieldConfig);

  // Data factors
  DATA_VECTOR(w_i);       	// weights for each obs
  DATA_VECTOR(l_i);       	// length for each obs
  DATA_FACTOR(X_i)          // Knot for each obs
  DATA_FACTOR(T_i)          // year for each obs
  DATA_FACTOR(S_i)          // sex for each obs
  DATA_FACTOR(D_i)          // date for each obs
  DATA_ARRAY(X_a_xjt);		    // Covariate design array (observation x covariate x year)
  DATA_ARRAY(X_b_xjt);		    // Covariate design array (observation x covariate x year)

  // Aniso objects
  DATA_INTEGER(n_tri);      //  Number of triangles
  DATA_VECTOR(Tri_Area);
  DATA_MATRIX(E0);
  DATA_MATRIX(E1);
  DATA_MATRIX(E2);
  DATA_FACTOR(TV);        //  This already includes the -1 for indexing in C
  DATA_SPARSE_MATRIX(G0_inv);

  // SPDE objects
  DATA_SPARSE_MATRIX(G0);
  DATA_SPARSE_MATRIX(G1);
  DATA_SPARSE_MATRIX(G2);

  // Fixed effects
  PARAMETER_VECTOR(ln_H_input); // Anisotropy parameters
  PARAMETER(alpha);   
  PARAMETER(beta);   
  PARAMETER(rho);             // Autocorrelation (i.e. density dependence)
  PARAMETER(gamma_a_sex);        // Covariate effect
  PARAMETER(gamma_b_sex);        // Covariate effect
  PARAMETER(gamma_a_date);        // Covariate effect
  PARAMETER(gamma_b_date);        // Covariate effect
  PARAMETER_VECTOR(gamma_a_j);        // Covariate effect
  PARAMETER_VECTOR(gamma_b_j);        // Covariate effect
  PARAMETER(log_sigma);
  PARAMETER(log_nu_E_a);      // nu = kappa*tau; margSD = 1 / (nu*sqrt(4*pi))
  PARAMETER(log_nu_O_a);      
  PARAMETER(log_SigmaD_a);
  PARAMETER(log_nu_E_b);      
  PARAMETER(log_nu_O_b);      
  PARAMETER(log_SigmaD_b);
  PARAMETER(log_kappa);      // Controls range of spatial variation

  // Random effects
  PARAMETER_VECTOR(Omega_input_a);   // Spatial variation in carrying capacity
  PARAMETER_ARRAY(Epsilon_input_a);  // Spatial process variation
  PARAMETER_VECTOR(Delta_a_t);
  PARAMETER_VECTOR(Omega_input_b);   // Spatial variation in carrying capacity
  PARAMETER_ARRAY(Epsilon_input_b);  // Spatial process variation
  PARAMETER_VECTOR(Delta_b_t);

  using namespace density;
  int i,j,s,t,x; // Observation, Covariate, Mesh-vertex, Year, Knot, 
  Type g = 0;
  
  // Derived values
  Type log_l_mean = sum( log(l_i) ) / l_i.size();
  Type kappa2 = exp(2.0*log_kappa);
  Type kappa4 = kappa2*kappa2;
  Type pi = 3.141592;
  Type Range = sqrt(8) / exp( log_kappa );
  Type SigmaM = exp(log_sigma);
  Type log_tau_E_a = log_nu_E_a - log_kappa;      // log-inverse SD of Epsilon
  Type log_tau_O_a = log_nu_O_a - log_kappa;      // log-inverse SD of Omega
  Type log_tau_E_b = log_nu_E_b - log_kappa;      // log-inverse SD of Epsilon
  Type log_tau_O_b = log_nu_O_b - log_kappa;      // log-inverse SD of Omega
  Type SigmaE_a = 1 / sqrt(4*pi*exp(2*log_tau_E_a)*exp(2*log_kappa));
  Type SigmaO_a = 1 / sqrt(4*pi*exp(2*log_tau_O_a)*exp(2*log_kappa));
  Type SigmaD_a = exp( log_SigmaD_a );
  Type SigmaE_b = 1 / sqrt(4*pi*exp(2*log_tau_E_b)*exp(2*log_kappa));
  Type SigmaO_b = 1 / sqrt(4*pi*exp(2*log_tau_O_b)*exp(2*log_kappa));
  Type SigmaD_b = exp( log_SigmaD_b );

  // Anisotropy elements
  matrix<Type> H(2,2);
  H(0,0) = exp(ln_H_input(0));
  H(1,0) = ln_H_input(1);
  H(0,1) = ln_H_input(1);
  H(1,1) = (1+ln_H_input(1)*ln_H_input(1)) / exp(ln_H_input(0));
  Type H_trace = H(0,0)+H(1,1);
  Type H_det = H(0,0)*H(1,1)-H(0,1)*H(1,0);
  Eigen::SparseMatrix<Type> G1_aniso(n_s,n_s); 
  Eigen::SparseMatrix<Type> G2_aniso(n_s,n_s); 
  // Calculate adjugate of H
  matrix<Type> adj_H(2,2);
  adj_H(0,0) = H(1,1);
  adj_H(0,1) = -1 * H(0,1);
  adj_H(1,0) = -1 * H(1,0);
  adj_H(1,1) = H(0,0);
  // Calculate new SPDE matrices
  if(Aniso==1){
    // Calculate G1 - pt. 1
    array<Type> Gtmp(n_tri,3,3);
    for(i=0; i<n_tri; i++){    
      // 1st line: E0(i,) %*% adjH %*% t(E0(i,)), etc.    
      Gtmp(i,0,0) = (E0(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E0(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
      Gtmp(i,0,1) = (E1(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E1(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));  
      Gtmp(i,0,2) = (E2(i,0)*(E0(i,0)*adj_H(0,0)+E0(i,1)*adj_H(1,0)) + E2(i,1)*(E0(i,0)*adj_H(0,1)+E0(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
      Gtmp(i,1,1) = (E1(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E1(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
      Gtmp(i,1,2) = (E2(i,0)*(E1(i,0)*adj_H(0,0)+E1(i,1)*adj_H(1,0)) + E2(i,1)*(E1(i,0)*adj_H(0,1)+E1(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
      Gtmp(i,2,2) = (E2(i,0)*(E2(i,0)*adj_H(0,0)+E2(i,1)*adj_H(1,0)) + E2(i,1)*(E2(i,0)*adj_H(0,1)+E2(i,1)*adj_H(1,1))) / (4*Tri_Area(i));
    }
    // Calculate G1 - pt. 2
    for(i=0; i<n_tri; i++){
      int i0 = i;
      int i1 = i + n_tri; 
      int i2 = i + 2*n_tri; 
      G1_aniso.coeffRef(TV(i1),TV(i0)) = G1_aniso.coeffRef(TV(i1),TV(i0)) + (Gtmp(i,0,1));  
      G1_aniso.coeffRef(TV(i0),TV(i1)) = G1_aniso.coeffRef(TV(i0),TV(i1)) + (Gtmp(i,0,1));  
      G1_aniso.coeffRef(TV(i2),TV(i1)) = G1_aniso.coeffRef(TV(i2),TV(i1)) + (Gtmp(i,1,2));  
      G1_aniso.coeffRef(TV(i1),TV(i2)) = G1_aniso.coeffRef(TV(i1),TV(i2)) + (Gtmp(i,1,2));  
      G1_aniso.coeffRef(TV(i2),TV(i0)) = G1_aniso.coeffRef(TV(i2),TV(i0)) + (Gtmp(i,0,2));  
      G1_aniso.coeffRef(TV(i0),TV(i2)) = G1_aniso.coeffRef(TV(i0),TV(i2)) + (Gtmp(i,0,2));  
      G1_aniso.coeffRef(TV(i0),TV(i0)) = G1_aniso.coeffRef(TV(i0),TV(i0)) + (Gtmp(i,0,0));  
      G1_aniso.coeffRef(TV(i1),TV(i1)) = G1_aniso.coeffRef(TV(i1),TV(i1)) + (Gtmp(i,1,1));  
      G1_aniso.coeffRef(TV(i2),TV(i2)) = G1_aniso.coeffRef(TV(i2),TV(i2)) + (Gtmp(i,2,2));  
    }
    G2_aniso = G1_aniso * G0_inv * G1_aniso; 
  }

  // Calculate fields
  matrix<Type> Epsilon_a_xt(n_x,n_t);
  vector<Type> Omega_a_x(n_x);
  matrix<Type> Epsilon_b_xt(n_x,n_t);
  vector<Type> Omega_b_x(n_x);
  matrix<Type> log_CF_xt(n_x,n_t);     // Random fields + covariates
  matrix<Type> log_CF_RF_xt(n_x,n_t);  // Just random fields
  matrix<Type> b_xt(n_x,n_t);          // Random fields + covariates
  matrix<Type> b_RF_xt(n_x,n_t);       // Just random fields
  for (int x=0;x<n_x;x++){
    Omega_a_x(x) = Omega_input_a(x) / exp(log_tau_O_a);
    Omega_b_x(x) = Omega_input_b(x) / exp(log_tau_O_b);
    for (int t=0;t<n_t;t++){ 
      Epsilon_a_xt(x,t) = Epsilon_input_a(x,t) / exp(log_tau_E_a);
      Epsilon_b_xt(x,t) = Epsilon_input_b(x,t) / exp(log_tau_E_b);
      // Condition factor
      log_CF_RF_xt(x,t) = alpha + Omega_a_x(x) + Epsilon_a_xt(x,t);     // /(1-rho)
      log_CF_xt(x,t) = log_CF_RF_xt(x,t) + Delta_a_t(t);
      for(int j=0;j<n_j_a;j++){
        log_CF_xt(x,t) += gamma_a_j(j)*X_a_xjt(x,j,t);
      }
      // Allometry
      b_RF_xt(x,t) = beta + Omega_b_x(x) + Epsilon_b_xt(x,t);     
      b_xt(x,t) = b_RF_xt(x,t) + Delta_b_t(t);
      for(int j=0;j<n_j_b;j++){
        b_xt(x,t) += gamma_b_j(j)*X_b_xjt(x,j,t);
      }
    }
  }
  
  // Likelihood contribution from fields
  Eigen::SparseMatrix<Type> Q;
  if(Aniso==0){
    Q = kappa4*G0 + Type(2.0)*kappa2*G1 + G2;
  }
  if(Aniso==1){
    Q = kappa4*G0 + Type(2.0)*kappa2*G1_aniso + G2_aniso;
  }
  //g += SEPARABLE(AR1(rho),GMRF(Q))(Epsilon_input);
  if(FieldConfig(0)==1) g += GMRF(Q)(Omega_input_a);
  if(FieldConfig(3)==1) g += GMRF(Q)(Omega_input_b);
  for(t=0;t<n_t;t++){
    if(FieldConfig(1)==1) g += GMRF(Q)(Epsilon_input_a.col(t));  
    if(FieldConfig(2)==1) g -= dnorm( Delta_a_t(t), Type(0.0), SigmaD_a, true );  
    if(FieldConfig(4)==1) g += GMRF(Q)(Epsilon_input_b.col(t));  
    if(FieldConfig(5)==1) g -= dnorm( Delta_b_t(t), Type(0.0), SigmaD_b, true );  
  }

  // Likelihood contribution from observations
  vector<Type> log_w_i_hat(n_i);
  for (int i=0;i<n_i;i++){                     //  + Epsilon_jt(J_i(i),T_i(i))/(1-rho)
    // Condition factor
    log_w_i_hat(i) = log_CF_xt(X_i(i),T_i(i));
    // Allometrry
    log_w_i_hat(i) += ( b_xt(X_i(i),T_i(i)) + gamma_b_sex*S_i(i) + gamma_b_date*D_i(i) ) * ( log(l_i(i)) - log_l_mean );
    // Sex and date
    log_w_i_hat(i) += gamma_a_sex*S_i(i) + gamma_a_date*D_i(i);
    g -= dnorm( log(w_i(i)), log_w_i_hat(i), SigmaM, true );
  }

  // Diagnostics
  REPORT( log_l_mean );
  // Spatial field summaries
  REPORT( H );
  REPORT( Range );
  REPORT( SigmaE_a );
  REPORT( SigmaO_a );
  REPORT( SigmaD_a );
  REPORT( SigmaE_b );
  REPORT( SigmaO_b );
  REPORT( SigmaD_b );
  REPORT( SigmaM );
  ADREPORT( Range );
  ADREPORT( SigmaE_a );
  ADREPORT( SigmaO_a );
  ADREPORT( SigmaD_a );
  ADREPORT( SigmaE_b );
  ADREPORT( SigmaO_b );
  ADREPORT( SigmaD_b );
  ADREPORT( SigmaM );
  // Fields
  REPORT( Epsilon_a_xt );
  REPORT( Omega_a_x );
  REPORT( Delta_a_t );
  REPORT( log_CF_xt );
  REPORT( log_CF_RF_xt );
  REPORT( Epsilon_b_xt );
  REPORT( Omega_b_x );
  REPORT( Delta_b_t );
  REPORT( b_xt );
  REPORT( b_RF_xt );
  return g;
}
