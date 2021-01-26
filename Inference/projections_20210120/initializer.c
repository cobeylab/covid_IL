// Number of sub-compartments per state 
int alpha_E_int = round(alpha_E);
int alpha_P_int = round(alpha_P); 
int alpha_IM_int = round(alpha_IM); 
int alpha_IS_int = round(alpha_IS); 
int alpha_IH1_int = round(alpha_IH1);
int alpha_IH2_int = round(alpha_IH2); 
int alpha_IC_int = round(alpha_IC);
int alpha_preIC_int = round(alpha_preIC);

// region-specific population size
const double *N = &N_1; //region-specific
const double *num_init = &num_init_1;

// initialize compartment pointers
double *S = &S_1_1;
double *R = &R_1_1;
double *D = &D_1_1;

double *E = &E1_1_1;
double *P = &P1_1_1;
double *IM = &IM1_1_1;
double *IM_dead = &IM_dead1_1_1;
double *IS = &IS1_1_1;

// Chains of severe infectious states
double *IH1 = &IH1_1_1_1;
double *IH2 = &IH2_1_1_1;
double *preIC = &preIC_1_1_1;
double *IC = &IC_1_1_1;

// Accumulators
double *new_deaths = &new_deaths_1_1;
double *new_hosp_deaths = &new_hosp_deaths_1_1;
double *new_nonhosp_deaths = &new_nonhosp_deaths_1_1;
double *new_hospitalizations = &new_hospitalizations_1_1;
double *Inc = &Inc_1_1;
double *DeathReportTrack = &DeathReportTrack_1_1;

// Draw initial age distribution using deaths before May 1 and IFR from Salje et al

for (int region=0; region<n_regions; region += 1){
      double excess = observed_all_cause - (exp_all_cause + 0);
      if (excess < 0){
          excess = 0;
      }
      double obs_nh = ll_deaths - ll_hosp;
      double exp_nh = excess - (ll_hosp  / reporting_cap);
      double nhd_report = obs_nh / exp_nh;
      if (nhd_report >= reporting_cap){            
        nhd_report = reporting_cap;
      } else if (nhd_report < 0 ){
        nhd_report = 0;  
      }
      DeathReportTrack[region] = nhd_report;

      double rN = rbinom(N[region], num_init[region]);
      int EStart = region * alpha_E_int;
      for (int j=0; j<alpha_E_int; j++){
        int offset = EStart + j;
        if (j == 0){
          E[offset] = rN;
        } 
        else{
          E[offset] = 0;    
        }
      }

      int PStart = region * alpha_P_int;
      for (int j=0; j<alpha_P_int; j++){
        int offset = PStart + j; 
        P[offset] = 0;
      }
      int IMStart = region * alpha_IM_int;
      for (int j=0; j<alpha_IM_int; j++){
        int offset = IMStart + j; 
        IM[offset] = 0;
        IM_dead[offset] = 0;
      }
      int ISStart = region * alpha_IS_int;
      for (int j=0; j<alpha_IS_int; j++){
        int offset = ISStart + j; 
        IS[offset] = 0;
      }
      int IH1Start = region * alpha_IH1_int;
      for (int j=0; j<alpha_IH1_int; j++){
        int offset = IH1Start + j; 
        IH1[offset] = 0;
      }
      int IH2Start = region * alpha_IH2_int;
      for (int j=0; j<alpha_IH2_int; j++){
        int offset = IH2Start + j; 
        IH2[offset] = 0;
      }
      int preIC_Start = region * alpha_preIC_int;
      for (int j=0; j<alpha_preIC_int; j++){
        int offset = preIC_Start + j;
        preIC[offset] = 0;
      }
      int IC_Start = region * alpha_IC_int;
      for (int j=0; j<alpha_IC_int; j++){
        int offset = IC_Start + j;
        IC[offset] = 0;
      }
      
      S[region] = N[region] - rN;
      R[region] = 0;
      D[region] = 0;
      new_deaths[region] = 0;
      new_hosp_deaths[region] = 0;
      new_nonhosp_deaths[region] = 0;
      new_hospitalizations[region] = 0;
      Inc[region] = rN;
}
