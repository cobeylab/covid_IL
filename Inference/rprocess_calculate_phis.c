
// Number of age groups 
int num_age_groups = n_age_groups;
//int dim_contact_matrix = n_age_groups*n_age_groups;

// Number of sub-compartments per state 
int alpha_E_int = alpha_E;
int alpha_A_int =  alpha_A;
int alpha_P_int = alpha_P; 
int alpha_IM_int = alpha_IM; 
int alpha_IS_int = alpha_IS; 
int alpha_IH1_int = alpha_IH1;
int alpha_IH2_int = alpha_IH2; 
int alpha_IH3_int = alpha_IH3;
int alpha_IH4_int = alpha_IH4;
int alpha_IC2_int = alpha_IC2; 
int alpha_IC3_int = alpha_IC3;


// add noise to beta for simulations
const int beta_noise = add_noise_to_beta;
const double amp_noise = beta_noise_amplitude;

// initialize parameter pointers
const double *rho = &rho_1;
const double *eta = &eta_1;
const double *kappa = &kappa_1;
const double *region_non_hosp = &region_non_hosp_1;
const double *zeta_s = &zeta_s_1;
const double *zeta_h = &zeta_h_1;
const double *zeta_c = &zeta_c_1;

const double *psi1 = &psi1_1;
const double *psi2 = &psi2_1;
const double *psi3 = &psi3_1;
const double *psi4 = &psi4_1;

const double *gamma_m = &gamma_m_1;
const double *gamma_h = &gamma_h_1;
const double *gamma_c = &gamma_c_1;

const double *mu_c = &mu_c_1;
const double *mu_h = &mu_h_1;
const double *mu_m = &mu_m_1;

const double *q = &q_1;
const double *age_beta_scales = &age_beta_scales_1;

// initialize contact matrix pointer for each type of contact, starts at the contact entry for 1, 1. Assumes that the contact rates are all in one continuous block of memory
const double *C_home = &C_home_1_1_1;    
const double *C_work = &C_work_1_1_1;  
const double *C_school = &C_school_1_1_1;  
const double *C_other = &C_other_1_1_1;  

// Initialize interventions pointer 
const double *scale_home_vec = &scale_home_1;
const double *scale_work_vec = &scale_work_1;
const double *scale_school_vec = &scale_school_1;
const double *scale_other_vec = &scale_other_1;

const double *scale_beta = &scale_beta_1;

const double *N = &N_1_1;
const double *beta2 = &beta2_1;
const double *beta1 = &beta1_1;
double *new_mild_infections = &new_mild_infections_1_1;
double *new_symptomatic_infections = &new_symptomatic_infections_1_1;
double *new_deaths = &new_deaths_1_1;
double *new_hosp_deaths = &new_hosp_deaths_1_1;
double *new_IH1 = &new_IH1_1_1;
double *new_IC2 = &new_IC2_1_1;
double *new_IC3 = &new_IC3_1_1;
double *new_IH4 = &new_IH4_1_1;
double *Inc = &Inc_1_1;
// initialize compartment pointers
double *S = &S_1_1;
double *R = &R_1_1;
double *D = &D_1_1;
double *E = &E1_1_1;
double *A = &A1_1_1;
double *P = &P1_1_1;
double *IM = &IM1_1_1;
double *IM_dead = &IM_dead1_1_1;
double *IS = &IS1_1_1;
// Chains of severe infectious states
double *IH1 = &IH1_1_1_1;
double *IH2 = &IH2_1_1_1;
double *IH3 = &IH3_1_1_1;
double *IH4 = &IH4_1_1_1;
double *IC2 = &IC2_1_1_1;
double *IC3 = &IC3_1_1_1;


// Code to have the option of looking at a single region
int start_loop;
int end_loop;
if (region_to_test < 0){
  start_loop = 0;
  end_loop = n_regions;
} else{
  start_loop = region_to_test - 1;
  end_loop = region_to_test;
}
for (int region=start_loop; region<end_loop; region += 1){
    double frac_of_deaths_non_hospitalized = region_non_hosp[region];

    for (int j=0; j<num_age_groups; j += 1){

      // double frac_of_deaths_non_hospitalized = runif(0.2, 0.4); // fraction of deaths that are non-hospitalized, draw it separately for every age group at every timestep

      // calculate lambda
      double lambda = 0;
      double betaT = beta1[region];

      // Nursing home scaling
      double C_nurse = age_beta_scales[j] * b_elderly;

      // Get the starting index (i.e., if you are on age group 2, you want to go to column 2 of the matrix)
      int cMatrixColStartIndex = j * num_age_groups + num_age_groups*num_age_groups*region;

      for (int i=0; i<num_age_groups; i+=1){
        int cMatrixOffset = cMatrixColStartIndex + i;
        double C = C_home[cMatrixOffset] + C_work[cMatrixOffset] + C_school[cMatrixOffset] + C_other[cMatrixOffset];

        // implement the intervention and sum contacts 
         double scale_home = scale_home_vec[region];
         double scale_work = scale_work_vec[region];
         double scale_school = scale_school_vec[region];
         double scale_other = scale_other_vec[region];
         
         C = C_home[cMatrixOffset]*scale_home + C_work[cMatrixOffset]*scale_work + C_school[cMatrixOffset]*scale_school + C_other[cMatrixOffset]*scale_other;




         if ((i == 7) | (i == 8) | (i == 9)){
            C += C_nurse;
         }

         
         if(use_post_intervention_beta > 0 ){
              betaT = beta2[region];
          }
          else{
            betaT = beta1[region];
          }
      
    
      betaT = betaT*scale_beta[region];
      
      // add noise to beta when evaluating increases in post-intervention transmission rate
      if(beta_noise == 1){
          betaT = rnorm(betaT, amp_noise*betaT);
          if(betaT < beta2[region]){
            betaT = beta2[region];
          }
      }
      
      betaT = betaT;
        
        double presymptomatic = 0; // This will be calculated by looping over all alpha_P subcompartments of P
        int PStart = i * alpha_P_int + region * alpha_P_int * num_age_groups;
        for (int x=0; x<alpha_P_int; x++){
          int offset = PStart + x;
          presymptomatic += P[offset];
        }
        
        double asymptomatic = 0; // This will be calculated by looping over all alpha_A subcompartments of A
        int AStart = i * alpha_A_int + region * alpha_A_int * num_age_groups;
        for (int x=0; x<alpha_A_int; x++){
          int offset = AStart + x;
          asymptomatic += A[offset];
        }
        
        double infectious_mild = 0; // This will be calculated by looping over all alpha_IM subcompartments of IM
        int IMStart = i * alpha_IM_int + region * alpha_IM_int * num_age_groups;
        for (int x=0; x<alpha_IM_int; x++){
          int offset = IMStart + x;
          infectious_mild += IM[offset] + IM_dead[offset];
        }
        
        double infectious_severe = 0; // This will be calculated by looping over all alpha_IS subcompartments of IS
        int ISStart = i * alpha_IS_int + region * alpha_IS_int * num_age_groups;
        for (int x=0; x<alpha_IS_int; x++){
          int offset = ISStart + x;
          infectious_severe += IS[offset];
        }

        lambda += C / N[i + (region * num_age_groups)] * betaT * q[j] *(presymptomatic + asymptomatic + infectious_mild + infectious_severe);
      }
      
      int S_index = j + (region * num_age_groups);
      // Outflow from S
      double dS = rbinom(S[S_index], 1-exp(-lambda * dt));
      S[S_index] += -dS;
      
      if(ISNA(S[S_index])){
        Rprintf("Region %d, S=%f, t=%f", region, S[S_index], t);
      }

      // Outflows from E
      int EStart = j * alpha_E_int + region * alpha_E_int * num_age_groups;
      double dE[alpha_E_int];
      double Pr_E = 1-exp(-sigma * alpha_E * dt);
      dE[0] = rbinom(E[EStart], Pr_E);
      E[EStart] += dS - dE[0];
      
      for (int i=1; i<alpha_E_int; i++){
        int offset = EStart + i;
        dE[i] = rbinom(E[offset], Pr_E);
        E[offset] += dE[i - 1] - dE[i];         
      }
      // Split the last outflow from E into Presymptomatic and Asymptomatic
      double dE_A = rbinom(dE[alpha_E_int-1], rho[j]);
      double dE_P = dE[alpha_E_int-1] - dE_A;
      
      
      // Outflows from A
      int AStart = j * alpha_A_int + region * alpha_A_int * num_age_groups;
      double dA[alpha_A_int];
      double Pr_A = 1-exp(-eta[j] * alpha_A * dt);
      dA[0] = rbinom(A[AStart], Pr_A);
      A[AStart] += dE_A - dA[0];
      for (int i=1; i<alpha_A_int; i++){
        int offset = AStart + i;
        dA[i] = rbinom(A[offset], Pr_A);
        A[offset] += dA[i - 1] - dA[i];
      }
      
      
      // Outflows from P
      int PStart = j * alpha_P_int + region * alpha_P_int * num_age_groups;
      double dP[alpha_P_int];
      double Pr_P = 1-exp(-zeta_s[j] * alpha_P * dt);
      dP[0] = rbinom(P[PStart], Pr_P);
      P[PStart] += dE_P - dP[0];
      for (int i=1; i<alpha_P_int; i++){
        int offset = PStart + i;
        dP[i] = rbinom(P[offset], Pr_P);
        P[offset] += dP[i-1] - dP[i];
      }
      // Split the last outflow from P into severe and mild
      double dP_s = rbinom(dP[alpha_P_int-1], kappa[j]);
      double dP_m = dP[alpha_P_int-1] - dP_s;
      
      //Split outflow from P_m into death and recovery
      double phi = 1 - (kappa[j] * (1 - psi1[j] - psi2[j]) * frac_of_deaths_non_hospitalized) / ((1 - frac_of_deaths_non_hospitalized) * (1 - kappa[j]));
      double dP_m_to_recover = rbinom(dP_m, phi); // phi is fraction mild cases that recover
      double dP_m_to_death = rbinom(dP_m, (1-phi));
      
      int IMStart = j * alpha_IM_int + region * alpha_IM_int * num_age_groups;
      double dIM[alpha_IM_int];
      double Pr_IM = 1-exp(-gamma_m[j] * alpha_IM * dt);
      dIM[0] = rbinom(IM[IMStart], Pr_IM);
      IM[IMStart] += dP_m_to_recover - dIM[0];
      for (int i=1; i<alpha_IM_int; i++){
        int offset = IMStart + i;
        dIM[i] = rbinom(IM[offset], Pr_IM);
        IM[offset] += dIM[i - 1] - dIM[i];
      }
      
      // Outflows from IM --> Death
      int IM_dead_Start = j * alpha_IM_int + region * alpha_IM_int * num_age_groups;
      double dIM_dead[alpha_IM_int];
      double Pr_IM_dead = 1-exp(-mu_m[j] * alpha_IM * dt); // mu_m is death rate for non-hospitalized deaths
      dIM_dead[0] = rbinom(IM_dead[IM_dead_Start], Pr_IM_dead);
      IM_dead[IM_dead_Start] += dP_m_to_death - dIM_dead[0];
      for (int i=1; i<alpha_IM_int; i++){
        int offset = IM_dead_Start + i;
        dIM_dead[i] = rbinom(IM_dead[offset], Pr_IM_dead);
        IM_dead[offset] += dIM_dead[i - 1] - dIM_dead[i];
      }
      
      // Severe, but unhospitalized
      int ISStart = j * alpha_IS_int + region * alpha_IS_int * num_age_groups;
      double dIS[alpha_IS_int];
      double Pr_IS = 1-exp(-zeta_h[j] * alpha_IS * dt);
      dIS[0] = rbinom(IS[ISStart], Pr_IS);
      IS[ISStart] += dP_s - dIS[0];
      for (int i=1; i<alpha_IS_int; i++){
        int offset = ISStart + i;
        dIS[i] = rbinom(IS[offset], Pr_IS);
        IS[offset] += dIS[i - 1] - dIS[i];
      }
      // Split last outflow of IS into: deaths, ICU and recover, and recover
      int Hosp_Split[4];
      double hosp_dist[4] = {psi1[j], psi2[j], psi3[j], psi4[j]};
      rmultinom(dIS[alpha_IS_int-1], hosp_dist, 4, Hosp_Split);
      double S_to_recover = Hosp_Split[0];
      double S_to_ICU_recover = Hosp_Split[1];
      double S_to_ICU_death = Hosp_Split[2];
      double S_to_hosp_death = Hosp_Split[3]; 
      //Rprintf("%f, %f, %f, %f\n", psi1[j], psi2[j], psi3[j], 1-psi1[j]-psi2[j]-psi3[j]);
      // Hospitalized -> recover
      // Outflows from IH1
      int IH1_Start = j * alpha_IH1_int + region * alpha_IH1_int * num_age_groups;
      double dIH1[alpha_IH1_int];
      double Pr_IH1 = 1-exp(-gamma_h[j] * alpha_IH1 * dt);
      dIH1[0] = rbinom(IH1[IH1_Start], Pr_IH1);
      IH1[IH1_Start] += S_to_recover - dIH1[0];
      for (int i=1; i<alpha_IH1_int; i++){
        int offset = IH1_Start + i;
        dIH1[i] = rbinom(IH1[offset], Pr_IH1);
        IH1[offset] += dIH1[i - 1] - dIH1[i];
      }
      
      // Hospitalized -> ICU -> Recover
      // Outflows from IH2
      int IH2_Start = j * alpha_IH2_int + region * alpha_IH2_int * num_age_groups;
      double dIH2[alpha_IH2_int];
      double Pr_IH2 = 1-exp(-zeta_c[j] * alpha_IH2 * dt);
      dIH2[0] = rbinom(IH2[IH2_Start], Pr_IH2);
      IH2[IH2_Start] += S_to_ICU_recover - dIH2[0];
      for (int i=1; i<alpha_IH2_int; i++){
        int offset = IH2_Start + i;
        dIH2[i] = rbinom(IH2[offset], Pr_IH2);
        IH2[offset] += dIH2[i - 1] - dIH2[i];
      }
      // Outflows from IC2
      int IC2_Start = j * alpha_IC2_int + region * alpha_IC2_int * num_age_groups;
      double dIC2[alpha_IC2_int];
      double Pr_IC2 = 1-exp(-gamma_c[j] * alpha_IC2 * dt);
      dIC2[0] = rbinom(IC2[IC2_Start], Pr_IC2);
      IC2[IC2_Start] += dIH2[alpha_IH2_int -1] - dIC2[0];
      for (int i=1; i<alpha_IC2_int; i++){
        int offset = IC2_Start + i;
        dIC2[i] = rbinom(IC2[offset], Pr_IC2);
        IC2[offset] += dIC2[i - 1] - dIC2[i];
      }
      
      // Hospitalized -> ICU -> Dead
      // Outflows from IH3
      int IH3_Start = j * alpha_IH3_int + region * alpha_IH3_int * num_age_groups;
      double dIH3[alpha_IH3_int];
      double Pr_IH3 = 1-exp(-zeta_c[j] * alpha_IH3 * dt);
      dIH3[0] = rbinom(IH3[IH3_Start], Pr_IH3);
      IH3[IH3_Start] += S_to_ICU_death - dIH3[0];
      for (int i=1; i<alpha_IH3_int; i++){
        int offset = IH3_Start + i;
        dIH3[i] = rbinom(IH3[offset], Pr_IH3);
        IH3[offset] += dIH3[i - 1] - dIH3[i];
      }
      // Outflows from IC3
      int IC3_Start = j * alpha_IC3_int + region * alpha_IC3_int * num_age_groups;
      double dIC3[alpha_IC3_int];
      double Pr_IC3 = 1-exp(-mu_c[j] * alpha_IC3 * dt);
      dIC3[0] = rbinom(IC3[IC3_Start], Pr_IC3);
      IC3[IC3_Start] += dIH3[alpha_IH3_int -1] - dIC3[0];
      for (int i=1; i<alpha_IC3_int; i++){
        int offset = IC3_Start + i;
        dIC3[i] = rbinom(IC3[offset], Pr_IC3);
        IC3[offset] += dIC3[i - 1] - dIC3[i];
      }
      
      
      // Hospitalized --> dead 
      int IH4_Start = j * alpha_IH4_int + region * alpha_IH4_int * num_age_groups;
      double dIH4[alpha_IH4_int];
      double Pr_IH4 = 1-exp(-mu_h[j] * alpha_IH4 * dt);
      dIH4[0] = rbinom(IH4[IH4_Start], Pr_IH4);
      IH4[IH4_Start] += S_to_hosp_death - dIH4[0];
      for (int i=1; i<alpha_IH4_int; i++){
        int offset = IH4_Start + i;
        dIH4[i] = rbinom(IH4[offset], Pr_IH4);
        IH4[offset] += dIH4[i - 1] - dIH4[i];
      }
      
      // Update recovered and dead
      R[j + (region * num_age_groups)] += dA[alpha_A_int-1] + dIM[alpha_IM_int-1] + dIH1[alpha_IH1_int-1] + dIC2[alpha_IC2_int-1];
      D[j + (region * num_age_groups)] += dIC3[alpha_IC3_int-1] + dIH4[alpha_IH4_int-1] + dIM_dead[alpha_IM_int-1];
      
      new_mild_infections[j + (region * num_age_groups)] += dP_m;
      new_symptomatic_infections[j + (region * num_age_groups)] += dP_s + dP_m;
      new_deaths[j + (region * num_age_groups)] += dIC3[alpha_IC3_int - 1] + dIH4[alpha_IH4_int - 1] + dIM_dead[alpha_IM_int-1];
      new_hosp_deaths[j + (region * num_age_groups)] += dIC3[alpha_IC3_int - 1] + dIH4[alpha_IH4_int - 1];
      Inc[j + (region * num_age_groups)] += dE[alpha_E_int - 1];
      new_IH1[j + (region * num_age_groups)] += S_to_recover;
      new_IC2[j + (region * num_age_groups)] += dIH2[alpha_IH2_int - 1];
      new_IC3[j + (region * num_age_groups)] += dIH3[alpha_IH3_int - 1];
      new_IH4[j + (region * num_age_groups)] += S_to_hosp_death;
    }
  }


