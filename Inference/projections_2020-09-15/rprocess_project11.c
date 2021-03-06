// Number of age groups 
int num_age_groups = round(n_age_groups);
int reg = round(region_to_test);

// Number of sub-compartments per state 
int alpha_E_int = round(alpha_E);
int alpha_A_int = round(alpha_A);
int alpha_P_int = round(alpha_P); 
int alpha_IM_int = round(alpha_IM); 
int alpha_IS_int = round(alpha_IS); 
int alpha_IH1_int = round(alpha_IH1);
int alpha_IH2_int = round(alpha_IH2); 
int alpha_IH3_int = round(alpha_IH3);
int alpha_IH4_int = round(alpha_IH4);
int alpha_IC2_int = round(alpha_IC2); 
int alpha_IC3_int = round(alpha_IC3);

// initialize pointers for age/region-specific parameters
const double *rho = &rho_1; // age-specific
const double *IHR_logit = &IHR_logit_1; // age-specific
const double *IHR_min = &IHR_min_1; // age-specific
const double *IHR_max = &IHR_max_1; // age-specific

//const double *kappa = &kappa_1; // age-specific
const double *q = &q_1; // age-specific
const double *age_beta_scales = &age_beta_scales_1; // age-specific

const double *region_non_hosp = &region_non_hosp_1; // region-specific
const double *beta2 = &beta2_1; //region-specific
const double *beta1 = &beta1_1; //Region-specific scaling on the maximum beta1

const double *N = &N_1_1; //age and region specific

// initialize contact matrix pointer for each type of contact, starts at the contact entry for 1, 1. Assumes that the contact rates are all in one continuous block of memory
  const double *C = &C_1_1_1;    

// Initialize interventions pointer 
  const double *t_phase3 = &t_phase3_1;
  const double *t_phase3_max = &t_phase3_max_1;
  const double *scale_phase3 = &scale_phase3_1;
  const double *t_phase4 = &t_phase4_1;
  const double *t_phase4_max = &t_phase4_max_1;
  const double *scale_phase4 = &scale_phase4_1;
  const double *t_phase5 = &t_phase5_1;
  const double *t_phase5_max = &t_phase5_max_1;
  const double *scale_phase5 = &scale_phase5_1;



// Age-specific hospitalization parameters

  const double *gamma_h = &gamma_h_1; // age-specific
  const double *zeta_c = &zeta_c_1; // age-specific
  const double *mu_c = &mu_c_1; // age-specific
  const double *psi1 = &psi1_1;
  const double *psi2 = &psi2_1;
  const double *psi3 = &psi3_1;
  const double *psi4 = &psi4_1;

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

// accumulators
double *new_symptomatic_infections = &new_symptomatic_infections_1_1;
double *new_deaths = &new_deaths_1_1;
double *new_hosp_deaths = &new_hosp_deaths_1_1;
double *new_nonhosp_deaths = &new_nonhosp_deaths_1_1;
double *new_hospitalizations = &new_hospitalizations_1_1;
double *Inc = &Inc_1_1;


// Code to have the option of looking at a single region
int start_loop;
int end_loop;
if (reg < 0){
  start_loop = 0;
  end_loop = n_regions;
} else{
  start_loop = reg - 1;
  end_loop = reg;
}
// Convert zeta_h logit to zeta_h
double zeta_h = (zeta_h_logit + zeta_h_min) * zeta_h_max / (1 + zeta_h_min);
int regparam = 0;

for (int region=start_loop; region<end_loop; region += 1){
  
    if (region == 0 | region == 1){
        regparam = 0;
    } else if (region == 2 | region == 5){
        regparam = 1;
    } else if (region == 3 | region == 4){
        regparam = 3;
    } else if (region == 10){
        regparam = 4;
    } else{
        regparam = 2;
    }

    // set region-specific params
      double frac_of_deaths_non_hospitalized = region_non_hosp[regparam];

    for (int j=0; j<num_age_groups; j += 1){

      // Calculate kappa from IHR
      double IHR = (IHR_logit[j] + IHR_min[j]) * IHR_max[j] / (1 + IHR_min[j]);
      double kappa = IHR / (1 - rho[j]);


      // calculate lambda
      double lambda = 0;

      // Get the starting index (i.e., if you are on age group 2, you want to go to column 2 of the matrix)
      int cMatrixColStartIndex = j * num_age_groups + num_age_groups*num_age_groups*regparam;

      for (int i=0; i<num_age_groups; i+=1){
        int cMatrixOffset = cMatrixColStartIndex + i;

         // Calculate transmission
           double transmission;
           double transmission_min = beta2[region];
           double transmission_max = beta1[region];

           int t_sip = 62;
           if (t < t_sip){ //pre-sip
             transmission = beta1[region];
           } else if (t < t_phase3[region]){ // phase 2 transmission
             double t_coord = t - t_sip + 1;
             double slope = (transmission_max - transmission_min) / 7;                 
             transmission = (t_coord < 7) ? transmission_max - slope * t_coord : transmission_min;
          } else if (t < t_phase4[region]){ // Phase 3 transmission
              double phase3_slope = (scale_phase3[region] * transmission_min) / (t_phase3_max[region] - t_phase3[region]);
              double transmission_phase3 = transmission_min + phase3_slope * (t - t_phase3[region]);
              double transmission_phase3_max = transmission_min + phase3_slope * (t_phase3_max[region] - t_phase3[region]);
              transmission = (t >= t_phase3_max[region]) ? transmission_phase3_max : transmission_phase3;
          } else if (t < t_phase5[region]){ // Phase 4 transmission
              double phase4_change = (scale_phase4[region] - scale_phase3[region]) * transmission_min;
              double phase4_starting_point = transmission_min * (1 + scale_phase3[region]);
              double phase4_slope = phase4_change / (t_phase4_max[region] - t_phase4[region]);
              double transmission_phase4 = phase4_starting_point + phase4_slope * (t - t_phase4[region]);
              double transmission_phase4_max = phase4_starting_point + phase4_slope * (t_phase4_max[region] - t_phase4[region]);
              transmission = (t >= t_phase4_max[region]) ? transmission_phase4_max : transmission_phase4;
          }else{
              double phase5_change = (scale_phase5[region] - scale_phase4[region]) * transmission_min;
              double phase5_starting_point = transmission_min * (1 + scale_phase4[region]);
              double phase5_slope = phase5_change / (t_phase5_max[region] - t_phase5[region]);
              double transmission_phase5 = phase5_starting_point + phase5_slope * (t - t_phase5[region]);
              double transmission_phase5_final = phase5_starting_point + phase5_slope * (t_phase5_max[region] - t_phase5[region]);
              transmission = (t >= t_phase5_max[region]) ? transmission_phase5_final : transmission_phase5; 
          }
            

        // Sum up everyone who's infectious
        double infectious = 0; // This will be calculated by looping over all alpha_P subcompartments of P
        int PStart = i * alpha_P_int + region * alpha_P_int * num_age_groups;
        for (int x=0; x<alpha_P_int; x++){
          int offset = PStart + x;
          infectious += P[offset];
        }

        int AStart = i * alpha_A_int + region * alpha_A_int * num_age_groups;
        for (int x=0; x<alpha_A_int; x++){
          int offset = AStart + x;
          infectious += A[offset];
        }

        int IMStart = i * alpha_IM_int + region * alpha_IM_int * num_age_groups;
        for (int x=0; x<alpha_IM_int; x++){
          int offset = IMStart + x;
          infectious += IM[offset] + IM_dead[offset];
        }

        int ISStart = i * alpha_IS_int + region * alpha_IS_int * num_age_groups;
        for (int x=0; x<alpha_IS_int; x++){
          int offset = ISStart + x;
          infectious += IS[offset];
        }

        lambda += transmission * C[cMatrixOffset] * q[j] * infectious / N[i + (region * num_age_groups)];
      }
      // Outflow from S
      double dS = rbinom(S[j + (region * num_age_groups)], 1-exp(-lambda * dt));
      S[j + (region * num_age_groups)] += -dS;

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
      double Pr_A = 1-exp(-eta * alpha_A * dt);
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
      double Pr_P = 1-exp(-zeta_s * alpha_P * dt);
      dP[0] = rbinom(P[PStart], Pr_P);
      P[PStart] += dE_P - dP[0];
      for (int i=1; i<alpha_P_int; i++){
        int offset = PStart + i;
        dP[i] = rbinom(P[offset], Pr_P);
        P[offset] += dP[i-1] - dP[i];
      }
      // Split the last outflow from P into severe and mild
      double dP_s = rbinom(dP[alpha_P_int-1], kappa);
      double dP_m = dP[alpha_P_int-1] - dP_s;
      
      //Split outflow from P_m into death and recovery
      double phi = 1 - (kappa * (1 - psi1[j] - psi2[j]) * frac_of_deaths_non_hospitalized) / ((1 - frac_of_deaths_non_hospitalized) * (1 - kappa));
      double dP_m_to_recover = rbinom(dP_m, phi); // phi is fraction mild cases that recover
      double dP_m_to_death = rbinom(dP_m, (1-phi));
      
      int IMStart = j * alpha_IM_int + region * alpha_IM_int * num_age_groups;
      double dIM[alpha_IM_int];
      double Pr_IM = 1-exp(-gamma_m * alpha_IM * dt);
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
      double Pr_IM_dead = 1-exp(-mu_m * alpha_IM * dt); // mu_m is death rate for non-hospitalized deaths
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
      double Pr_IS = 1-exp(-zeta_h * alpha_IS * dt);
      dIS[0] = rbinom(IS[ISStart], Pr_IS);
      IS[ISStart] += dP_s - dIS[0];
      for (int i=1; i<alpha_IS_int; i++){
        int offset = ISStart + i;
        dIS[i] = rbinom(IS[offset], Pr_IS);
        IS[offset] += dIS[i - 1] - dIS[i];
      }

      // Split last outflow of IS into: deaths, ICU and recover, and recover
      int Hosp_Split[4];
      double psisum = psi1[j] + psi2[j] + psi3[j] + psi4[j];
      double hosp_dist[4] = {psi1[j]/psisum, psi2[j]/psisum, psi3[j]/psisum, psi4[j]/psisum};
      rmultinom(dIS[alpha_IS_int-1], hosp_dist, 4, Hosp_Split);
      double S_to_recover = Hosp_Split[0];
      double S_to_ICU_recover = Hosp_Split[1];
      double S_to_ICU_death = Hosp_Split[2];
      double S_to_hosp_death = Hosp_Split[3]; 

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
      double Pr_IC2 = 1-exp(-gamma_c * alpha_IC2 * dt);
      dIC2[0] = rbinom(IC2[IC2_Start], Pr_IC2);
      IC2[IC2_Start] += dIH2[alpha_IH2_int -1] - dIC2[0];
      for (int i=1; i<alpha_IC2_int; i++){
        int offset = IC2_Start + i;
        dIC2[i] = rbinom(IC2[offset], Pr_IC2);
        IC2[offset] += dIC2[i - 1] - dIC2[i];
      }

      // Add recovered ICU patients back to the hospital
      IH1[IH1_Start] += dIC2[alpha_IC2_int - 1];
      
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
      double Pr_IH4 = 1-exp(-mu_h * alpha_IH4 * dt);
      dIH4[0] = rbinom(IH4[IH4_Start], Pr_IH4);
      IH4[IH4_Start] += S_to_hosp_death - dIH4[0];
      for (int i=1; i<alpha_IH4_int; i++){
        int offset = IH4_Start + i;
        dIH4[i] = rbinom(IH4[offset], Pr_IH4);
        IH4[offset] += dIH4[i - 1] - dIH4[i];
      }
      
      // Update recovered and dead
      R[j + (region * num_age_groups)] += dA[alpha_A_int-1] + dIM[alpha_IM_int-1] + dIH1[alpha_IH1_int-1];
      D[j + (region * num_age_groups)] += dIC3[alpha_IC3_int-1] + dIH4[alpha_IH4_int-1] + dIM_dead[alpha_IM_int-1];
      

      new_symptomatic_infections[j + (region * num_age_groups)] += dP_s + dP_m;

      new_deaths[j + (region * num_age_groups)] += dIC3[alpha_IC3_int - 1] + dIH4[alpha_IH4_int - 1] + dIM_dead[alpha_IM_int-1];
      new_hosp_deaths[j + (region * num_age_groups)] += dIC3[alpha_IC3_int - 1] + dIH4[alpha_IH4_int - 1];
      new_nonhosp_deaths[j + (region * num_age_groups)] += dIM_dead[alpha_IM_int-1];
      new_hospitalizations[j + (region * num_age_groups)] += dIS[alpha_IS_int-1];
      Inc[j + (region * num_age_groups)] += dE[alpha_E_int - 1];

      //new_mild_infections[j + (region * num_age_groups)] += dP_m;
      //new_IH1[j + (region * num_age_groups)] += S_to_recover;
      //new_IC2[j + (region * num_age_groups)] += dIH2[alpha_IH2_int - 1];
      //new_IC3[j + (region * num_age_groups)] += dIH3[alpha_IH3_int - 1];
      //new_IH4[j + (region * num_age_groups)] += S_to_hosp_death;
    }
  }


