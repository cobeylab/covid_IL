
// Number of age groups 
int num_age_groups = n_age_groups;
//int dim_contact_matrix = n_age_groups*n_age_groups;

// Number of interventions
int num_interventions = n_interventions;

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



// initialize parameter pointers
const double *rho = &rho_1;
const double *eta = &eta_1;
const double *kappa = &kappa_1;

const double *zeta_s = &zeta_s_1;
const double *zeta_h = &zeta_h_1;
const double *zeta_c = &zeta_c_1;

const double *psi1 = &psi1_1;
const double *psi2 = &psi2_1;
const double *psi3 = &psi3_1;

const double *gamma_m = &gamma_m_1;
const double *gamma_h = &gamma_h_1;
const double *gamma_c = &gamma_c_1;

const double *mu_c = &mu_c_1;
const double *mu_h = &mu_h_1;

// initialize contact matrix pointer for each type of contact, starts at the contact entry for 1, 1. Assumes that the contact rates are all in one continuous block of memory
const double *C_home = &C_home_1_1;    
const double *C_work = &C_work_1_1;  
const double *C_school = &C_school_1_1;  
const double *C_other = &C_other_1_1;  

// Initialize interventions pointer 
const double *t_start_intervention = &t_start_intervention_1;
const double *t_end_intervention = &t_end_intervention_1;
const double *scale_work_intervention = &scale_work_intervention_1;
const double *scale_home_intervention =  &scale_home_intervention_1;
const double *scale_school_intervention = &scale_school_intervention_1;
const double *scale_other_intervention = &scale_other_intervention_1;

const double *N = &N_1_1;
//const double *beta2 = &beta2_1;
double *new_mild_infections = &new_mild_infections_1_1;
double *new_severe_infections = &new_severe_infections_1_1;
double *new_deaths = &new_deaths_1_1;
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
double *IS = &IS1_1_1;
// Chains of severe infectious states
double *IH1 = &IH1_1_1_1;
double *IH2 = &IH2_1_1_1;
double *IH3 = &IH3_1_1_1;
double *IH4 = &IH4_1_1_1;
double *IC2 = &IC2_1_1_1;
double *IC3 = &IC3_1_1_1;


// loop over every age group
for (int region=0; region<n_regions; region +=1)
{
    double beta2;
    if (region == 0){
      beta2 = beta2_1;
    } else if (region ==1){
      beta2 = beta2_2;
    } else if (region ==2){
      beta2 = beta2_3;
    }
    for (int j=0; j<num_age_groups; j += 1){


      // calculate lambda
      double lambda = 0;
      double betaT = beta1;

      // Get the starting index (i.e., if you are on age group 2, you want to go to column 2 of the matrix)
      int cMatrixColStartIndex = j * num_age_groups;

      for (int i=0; i<num_age_groups; i+=1){
        int cMatrixOffset = cMatrixColStartIndex + i;
        double C = C_home[cMatrixOffset] + C_work[cMatrixOffset] + C_school[cMatrixOffset] + C_other[cMatrixOffset];

        // implement the intervention and sum contacts 
        for( int z = 0; z < num_interventions; z ++){
         double t_start = t_start_intervention[z];
         double t_end = t_end_intervention[z];
         double scale_home = scale_home_intervention[z];
         double scale_work = scale_work_intervention[z];
         double scale_school = scale_school_intervention[z];
         double scale_other = scale_other_intervention[z];
         
         if((t >= t_start) & (t<t_end)){
              betaT = beta2;
              C = C_home[cMatrixOffset]*scale_home + C_work[cMatrixOffset]*scale_work + C_school[cMatrixOffset]*scale_school + C_other[cMatrixOffset]*scale_other;
             }
          else{
            betaT = beta1;
          }
      }

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
          infectious_mild += IM[offset];
        }
        
        double infectious_severe = 0; // This will be calculated by looping over all alpha_IS subcompartments of IS
        int ISStart = i * alpha_IS_int + region * alpha_IS_int * num_age_groups;
        for (int x=0; x<alpha_IS_int; x++){
          int offset = ISStart + x;
          infectious_severe += IS[offset];
        }

        lambda += C / N[i + (region * num_age_groups)] * betaT * (presymptomatic + asymptomatic + infectious_mild + infectious_severe);
      }
      
      int S_index = j + (region * num_age_groups);
      // Outflow from S
      double dS = rbinom(S[S_index], 1-exp(-lambda * dt));
      S[S_index] += -dS;
      
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
      
      
      // Outflows from IM, assume no one dies.
      int IMStart = j * alpha_IM_int + region * alpha_IM_int * num_age_groups;
      double dIM[alpha_IM_int];
      double Pr_IM = 1-exp(-gamma_m[j] * alpha_IM * dt);
      dIM[0] = rbinom(IM[IMStart], Pr_IM);
      IM[IMStart] += dP_m - dIM[0];
      for (int i=1; i<alpha_IM_int; i++){
        int offset = IMStart + i;
        dIM[i] = rbinom(IM[offset], Pr_IM);
        IM[offset] += dIM[i - 1] - dIM[i];
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
      double hosp_dist[4] = {psi1[j], psi2[j], psi3[j], 1-psi1[j]-psi2[j]-psi3[j]};
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
      D[j + (region * num_age_groups)] += dIC3[alpha_IC3_int-1] + dIH4[alpha_IH4_int-1];

      new_mild_infections[j + (region * num_age_groups)] += dP_m;
      new_severe_infections[j + (region * num_age_groups)] += S_to_recover + S_to_ICU_recover;
      new_deaths[j + (region * num_age_groups)] += dIC3[alpha_IC3_int - 1] + dIH4[alpha_IH4_int - 1];
      Inc[j + (region * num_age_groups)] += dE[alpha_E_int - 1];
      new_IH1[j + (region * num_age_groups)] += S_to_recover;
      new_IC2[j + (region * num_age_groups)] += dIH2[alpha_IH2_int - 1];
      new_IC3[j + (region * num_age_groups)] += dIH3[alpha_IH3_int - 1];
      new_IH4[j + (region * num_age_groups)] += S_to_hosp_death;

    }
}