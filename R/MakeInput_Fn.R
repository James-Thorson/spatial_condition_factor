MakeInput_Fn <-
function( n_x, Version, MeshList, Data_Geostat, FieldConfig, ObsConfig, ObsModel, Aniso, NN, X_a_xjt, X_b_xjt ){

  if(Version%in%c("growth_v1a","growth_v1b","growth_v1c","growth_v2a","growth_v2a_fix","growth_v2a_bridge","growth_v2b","growth_v2c","growth_v2d","growth_v2e")){
    # Data
    if( Version %in% c("growth_v1a") ) TmbData = list("n_i"=nrow(Data_Geostat), "n_s"=MeshList$spde$n.spde, "n_j"=n_x, , "n_t"=diff(range(unique(Data_Geostat[,'Year'])))+1, "Aniso"=Aniso, "ObsModel"=ObsModel, "w_i"=Data_Geostat[,'WEIGHT_KG'], "l_i"=Data_Geostat[,'LENGTH_CM'], "J_i"=NN$nn.idx[,1]-1, "T_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=as(diag(1/diag(MeshList$spde$param.inla$M0)),"dgTMatrix"), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if( Version %in% c("growth_v1b") ) TmbData = list("n_i"=nrow(Data_Geostat), "n_s"=MeshList$spde$n.spde, "n_j"=n_x, "n_t"=diff(range(unique(Data_Geostat[,'Year'])))+1, "Aniso"=Aniso, "ObsModel"=ObsModel, "w_i"=Data_Geostat[,'WEIGHT_KG'], "l_i"=Data_Geostat[,'LENGTH_CM'], "J_i"=NN$nn.idx[,1]-1, "T_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "S_i"=sapply(Data_Geostat[,'Sex'],FUN=switch,0,1,0.5), "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=as(diag(1/diag(MeshList$spde$param.inla$M0)),"dgTMatrix"), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if( Version %in% c("growth_v1c","growth_v2a") ) TmbData = list("n_i"=nrow(Data_Geostat), "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=diff(range(unique(Data_Geostat[,'Year'])))+1, "n_j"=dim(X_a_xjt)[2], "Aniso"=Aniso, "ObsModel"=ObsModel, "w_i"=Data_Geostat[,'WEIGHT_KG'], "l_i"=Data_Geostat[,'LENGTH_CM'], "X_i"=NN$nn.idx[,1]-1, "T_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "S_i"=sapply(Data_Geostat[,'Sex'],FUN=switch,0,1,0.5), "X_xjt"=X_a_xjt, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=as(diag(1/diag(MeshList$spde$param.inla$M0)),"dgTMatrix"), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if( Version %in% c("growth_v2a_fix") ) TmbData = list("n_i"=nrow(Data_Geostat), "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=diff(range(unique(Data_Geostat[,'Year'])))+1, "n_j"=dim(X_a_xjt)[2], "Aniso"=Aniso, "ObsModel"=ObsModel, "FieldConfig"=FieldConfig, "w_i"=Data_Geostat[,'WEIGHT_KG'], "l_i"=Data_Geostat[,'LENGTH_CM'], "X_i"=NN$nn.idx[,1]-1, "T_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "S_i"=sapply(Data_Geostat[,'Sex'],FUN=switch,0,1,0.5), "X_xjt"=X_a_xjt, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=as(diag(1/diag(MeshList$spde$param.inla$M0)),"dgTMatrix"), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if( Version %in% c("growth_v2a_bridge") ) TmbData = list("n_i"=nrow(Data_Geostat), "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=diff(range(unique(Data_Geostat[,'Year'])))+1, "n_j"=dim(X_a_xjt)[2], "Aniso"=Aniso, "ObsModel"=ObsModel, "FieldConfig"=FieldConfig, "w_i"=Data_Geostat[,'WEIGHT_KG'], "l_i"=Data_Geostat[,'LENGTH_CM'], "X_i"=NN$nn.idx[,1]-1, "T_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "S_i"=sapply(Data_Geostat[,'Sex'],FUN=switch,0,1,0.5), "X_xjt"=X_a_xjt, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=as(diag(1/diag(MeshList$spde$param.inla$M0)),"dgTMatrix"), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if( Version %in% c("growth_v2b") ) TmbData = list("n_i"=nrow(Data_Geostat), "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=diff(range(unique(Data_Geostat[,'Year'])))+1, "n_j_a"=dim(X_a_xjt)[2], "n_j_b"=dim(X_b_xjt)[2], "Aniso"=Aniso, "ObsModel"=ObsModel, "FieldConfig"=FieldConfig, "w_i"=Data_Geostat[,'WEIGHT_KG'], "l_i"=Data_Geostat[,'LENGTH_CM'], "X_i"=NN$nn.idx[,1]-1, "T_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "S_i"=sapply(Data_Geostat[,'Sex'],FUN=switch,0,1,0.5), "X_a_xjt"=X_a_xjt, "X_b_xjt"=X_b_xjt, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=as(diag(1/diag(MeshList$spde$param.inla$M0)),"dgTMatrix"), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    if( Version %in% c("growth_v2c","growth_v2d","growth_v2e") ) TmbData = list("n_i"=nrow(Data_Geostat), "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=diff(range(unique(Data_Geostat[,'Year'])))+1, "n_j_a"=dim(X_a_xjt)[2], "n_j_b"=dim(X_b_xjt)[2], "Aniso"=Aniso, "ObsModel"=ObsModel, "FieldConfig"=FieldConfig, "w_i"=Data_Geostat[,'WEIGHT_KG'], "l_i"=Data_Geostat[,'LENGTH_CM'], "X_i"=NN$nn.idx[,1]-1, "T_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "S_i"=as.numeric(sapply(as.character(Data_Geostat[,'Sex']),FUN=switch,"f"=0,"m"=1,"u"=0.5)), "D_i"=(Data_Geostat[,'Date']-mean(Data_Geostat[,'Date']))/sd(Data_Geostat[,'Date']), "X_a_xjt"=X_a_xjt, "X_b_xjt"=X_b_xjt, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=as(diag(1/diag(MeshList$spde$param.inla$M0)),"dgTMatrix"), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2)
    # Parameters
    if( Version %in% c("growth_v1a") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=0, "beta"=3, "rho"=0, "log_sigma"=5, "log_tau_E"=log(1), "log_tau_O"=log(1), "log_kappa"=log(1), "Omega_input"=rep(0,TmbData$n_s), "Epsilon_input"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t))
    if( Version %in% c("growth_v1b") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=0, "beta"=3, "rho"=0, "gamma_sex"=0, "log_sigma"=5, "log_tau_E"=log(1), "log_tau_O"=log(1), "log_kappa"=log(1), "Omega_input"=rep(0,TmbData$n_s), "Epsilon_input"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t))
    if( Version %in% c("growth_v1c") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=0, "beta"=3, "rho"=0, "gamma_sex"=0, "gamma_j"=rep(0,TmbData$n_j), "log_sigma"=5, "log_tau_E"=log(1), "log_tau_O"=log(1), "log_kappa"=log(1), "Omega_input"=rep(0,TmbData$n_s), "Epsilon_input"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t))
    if( Version %in% c("growth_v2a","growth_v2a_fix") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=0, "beta"=3, "rho"=0, "gamma_sex"=0, "gamma_j"=rep(0,TmbData$n_j), "log_sigma"=5, "log_tau_E_a"=log(1), "log_tau_O_a"=log(1), "log_tau_E_b"=log(1), "log_tau_O_b"=log(1), "log_kappa"=log(1), "Omega_input_a"=rep(0,TmbData$n_s), "Epsilon_input_a"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t), "Omega_input_b"=rep(0,TmbData$n_s), "Epsilon_input_b"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t))
    if( Version %in% c("growth_v2a_bridge") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=mean(log(TmbData$l_i))*3, "beta"=3, "rho"=0, "gamma_sex"=0, "gamma_j"=rep(0,TmbData$n_j), "log_sigma"=5, "log_nu_E_a"=log(1), "log_nu_O_a"=log(1), "log_nu_E_b"=log(1), "log_nu_O_b"=log(1), "log_kappa"=log(1), "Omega_input_a"=rep(0,TmbData$n_s), "Epsilon_input_a"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t), "Omega_input_b"=rep(0,TmbData$n_s), "Epsilon_input_b"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t))
    if( Version %in% c("growth_v2b") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=0, "beta"=3, "rho"=0, "gamma_sex"=0, "gamma_a_j"=rep(0,TmbData$n_j_a), "gamma_b_j"=rep(0,TmbData$n_j_b), "log_sigma"=5, "log_nu_E_a"=log(1), "log_nu_O_a"=log(1), "log_nu_E_b"=log(1), "log_nu_O_b"=log(1), "log_kappa"=log(1), "Omega_input_a"=rep(0,TmbData$n_s), "Epsilon_input_a"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t), "Omega_input_b"=rep(0,TmbData$n_s), "Epsilon_input_b"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t))
    if( Version %in% c("growth_v2c","growth_v2d") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=mean(log(TmbData$l_i))*3, "beta"=3, "rho"=0, "gamma_a_sex"=0, "gamma_b_sex"=0, "gamma_a_date"=0, "gamma_b_date"=0, "gamma_a_j"=rep(0,TmbData$n_j_a), "gamma_b_j"=rep(0,TmbData$n_j_b), "log_sigma"=5, "log_nu_E_a"=log(1), "log_nu_O_a"=log(1), "log_nu_E_b"=log(1), "log_nu_O_b"=log(1), "log_kappa"=log(1), "Omega_input_a"=rep(0,TmbData$n_s), "Epsilon_input_a"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t), "Omega_input_b"=rep(0,TmbData$n_s), "Epsilon_input_b"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t))
    if( Version %in% c("growth_v2e") ) Parameters = list("ln_H_input"=c(0,0), "alpha"=mean(log(TmbData$l_i))*3, "beta"=3, "rho"=0, "gamma_a_sex"=0, "gamma_b_sex"=0, "gamma_a_date"=0, "gamma_b_date"=0, "gamma_a_j"=rep(0,TmbData$n_j_a), "gamma_b_j"=rep(0,TmbData$n_j_b), "log_sigma"=5, "log_nu_E_a"=log(1), "log_nu_O_a"=log(1), "log_SigmaD_a"=log(1), "log_nu_E_b"=log(1), "log_nu_O_b"=log(1), "log_SigmaD_b"=log(1), "log_kappa"=log(1), "Omega_input_a"=rep(0,TmbData$n_s), "Epsilon_input_a"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t), "Delta_a_t"=rep(0,TmbData$n_t), "Omega_input_b"=rep(0,TmbData$n_s), "Epsilon_input_b"=matrix(0,nrow=TmbData$n_s,ncol=TmbData$n_t), "Delta_b_t"=rep(0,TmbData$n_t))
    # Which are random
    if(Version %in% c("growth_v1a","growth_v1b","growth_v1c")) Random = c("Omega_input", "Epsilon_input")
    if(Version %in% c("growth_v2a","growth_v2a_fix","growth_v2a_bridge","growth_v2b","growth_v2c","growth_v2d")) Random = c("Omega_input_a", "Epsilon_input_a", "Omega_input_b", "Epsilon_input_b")
    if(Version %in% c("growth_v2e")) Random = c("Omega_input_a", "Epsilon_input_a", "Delta_a_t", "Omega_input_b", "Epsilon_input_b", "Delta_b_t")
    # Turn off parameters
    Map = list()
    Map[["rho"]] = factor(NA)
    if(Version %in% c("growth_v1c") ){
      if(FieldConfig['Omega_a']==0){
        Map[["Omega_input"]] = factor(rep(NA,TmbData$n_s))
        Map[["log_tau_O"]] = factor(NA)
      }
      if(FieldConfig['Epsilon_a']==0){
        Map[["Epsilon_input"]] = factor(matrix(NA,nrow=TmbData$n_s,ncol=TmbData$n_t))
        Map[["log_tau_E"]] = factor(NA)
      }
      if(FieldConfig['Omega_a']==0 & FieldConfig['Epsilon_a']==0) Map[["log_kappa"]] = factor(NA)
    }
    if(Version %in% c("growth_v2a","growth_v2a_fix","growth_v2a_bridge") ){
      if( sum(CovConfig_a)==0 & CovConception==FALSE ){
        Map[["gamma_j"]] = factor(NA)
      }
    }
    if(Version %in% c("growth_v2a","growth_v2a_fix","growth_v2a_bridge","growth_v2b") ){
      if( ObsConfig['Sex_a']==0){
        Map[["gamma_sex"]] = factor(NA)
      }
    }
    if(Version %in% c("growth_v2a","growth_v2a_fix","growth_v2a_bridge","growth_v2b","growth_v2c","growth_v2d","growth_v2e") ){
      if(FieldConfig['Omega_a']==0){
        Map[["Omega_input_a"]] = factor(rep(NA,TmbData$n_s))
        if(Version %in% c("growth_v2a","growth_v2a_fix") ) Map[["log_tau_O_a"]] = factor(NA)
        if(Version %in% c("growth_v2a_bridge","growth_v2b","growth_v2c","growth_v2d","growth_v2e") ) Map[["log_nu_O_a"]] = factor(NA)
      }
      if(FieldConfig['Epsilon_a']==0){
        Map[["Epsilon_input_a"]] = factor(matrix(NA,nrow=TmbData$n_s,ncol=TmbData$n_t))
        if(Version %in% c("growth_v2a","growth_v2a_fix") ) Map[["log_tau_E_a"]] = factor(NA)
        if(Version %in% c("growth_v2a_bridge","growth_v2b","growth_v2c","growth_v2d","growth_v2e") ) Map[["log_nu_E_a"]] = factor(NA)
      }
      if(FieldConfig['Omega_b']==0){
        Map[["Omega_input_b"]] = factor(rep(NA,TmbData$n_s))
        if(Version %in% c("growth_v2a","growth_v2a_fix") ) Map[["log_tau_O_b"]] = factor(NA)
        if(Version %in% c("growth_v2a_bridge","growth_v2b","growth_v2c","growth_v2d","growth_v2e") ) Map[["log_nu_O_b"]] = factor(NA)
      }
      if(FieldConfig['Epsilon_b']==0){
        Map[["Epsilon_input_b"]] = factor(matrix(NA,nrow=TmbData$n_s,ncol=TmbData$n_t))
        if(Version %in% c("growth_v2a","growth_v2a_fix") ) Map[["log_tau_E_b"]] = factor(NA)
        if(Version %in% c("growth_v2a_bridge","growth_v2b","growth_v2c","growth_v2d","growth_v2e") ) Map[["log_nu_E_b"]] = factor(NA)
      }
      if( all(FieldConfig[c('Omega_a','Epsilon_a','Omega_b','Epsilon_b')]==0) ) Map[["log_kappa"]] = factor(NA)
    }
    if(Version %in% c("growth_v2b","growth_v2c","growth_v2d","growth_v2e") ){
      if( sum(CovConfig_a)==0 & CovConception==FALSE ){
        Map[["gamma_a_j"]] = factor(NA)
      }
      if( sum(CovConfig_b)==0 & CovConception==FALSE ){
        Map[["gamma_b_j"]] = factor(NA)
      }
    }
    if(Version %in% c("growth_v2c","growth_v2d","growth_v2e") ){
      if( ObsConfig['Sex_a']==0 ) Map[["gamma_a_sex"]] = factor(NA)
      if( ObsConfig['Sex_b']==0 ) Map[["gamma_b_sex"]] = factor(NA)
      if( ObsConfig['Date_a']==0 ) Map[["gamma_a_date"]] = factor(NA)
      if( ObsConfig['Date_b']==0 ) Map[["gamma_b_date"]] = factor(NA)
    }
    if(Version %in% c("growth_v2e") ){
      if( FieldConfig['Delta_a']==0 ){
        Map[["Delta_a_t"]] = factor( rep(NA,length(Parameters[['Delta_a_t']])) )
        Map[["log_SigmaD_a"]] = factor( NA )
      }
      if( FieldConfig['Delta_b']==0 ){
        Map[["Delta_b_t"]] = factor( rep(NA,length(Parameters[['Delta_b_t']])) )
        Map[["log_SigmaD_b"]] = factor( NA )
      }
    }
    if(Aniso==0 | all(FieldConfig==0)) Map[['ln_H_input']] = factor( rep(NA,2) )
  }

  Return = list( "TmbData"=TmbData, "Parameters"=Parameters, "Random"=Random, "Map"=Map)
  return( Return )

}
