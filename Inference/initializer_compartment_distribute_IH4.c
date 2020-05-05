
// number of age groups, passed as parameter
int num_age_groups = n_age_groups;
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

//const double *num_init = &num_init_1;

int num_init_compartments = 3;


// initialize parameter pointers
const double *N = &N_1_1;
const double *age_dist = &age_dist_1_1;

const double *rho = &rho_1;
const double *kappa = &kappa_1;


// For initializing new cases;
double *new_mild_infections = &new_mild_infections_1_1;
double *new_symptomatic_infections = &new_symptomatic_infections_1_1;
double *new_deaths = &new_deaths_1_1;
double *new_IH1 = &new_IH1_1_1;
double *new_IC2 = &new_IC2_1_1;
double *new_IC3 = &new_IC3_1_1;
double *new_IH4 = &new_IH4_1_1;


// Right now just initialized at whatever values you want.
// Initial conditions (e.g. S0) are vectorized parameters.


for (int region=0; region<n_regions; region += 1 ){
  double num_init;

  if (region == 0){
    num_init = num_init_1;
  } else if (region==1){
    num_init = num_init_2;
  } else if (region==2){
    num_init = num_init_3;
  }

  // initialize multinomial draw to determine age distribution
  int rN[num_age_groups];
  double age_dist_region[num_age_groups];
  for (int ag=0; ag<num_age_groups; ag+=1){
    age_dist_region[ag] = age_dist[ag + region * num_age_groups];
  }
  rmultinom(ceil(num_init), age_dist_region, num_age_groups, rN);

  for (int i=0; i<num_age_groups; i += 1){
      double new_mild = 0;
      double new_asymp = 0;
      double new_IS = 0;

      // do the multinomial draw to determine compartments for age group i
      int compartment_init[num_init_compartments];
      double compartment_dist[num_init_compartments];
      
      compartment_dist[0] = rho[i]; // asymptomatic 
      compartment_dist[1] = (1-rho[i])*(1-kappa[i]); // mild infections
      compartment_dist[2] = (1-rho[i])*kappa[i]; // severe infections
      
      rmultinom(rN[i], compartment_dist, num_init_compartments, compartment_init);
     // Rprintf(\"asymp is %f\\n\",compartment_dist[0]);

      int EStart = i * alpha_E_int + region * alpha_E_int * num_age_groups;
      for (int j=0; j<alpha_E_int; j++){
        int offset = EStart + j; 
        E[offset] = 0;
      }

      int AStart = i * alpha_A_int + region * alpha_A_int * num_age_groups;
      for (int j=0; j<alpha_A_int; j++){
        int offset = AStart + j; 
        if (j == 0){
           A[offset] = compartment_init[0];
          }
         else{
          A[offset] = 0;
        }
         new_asymp += A[offset];
      }

      int PStart = i * alpha_P_int + region * alpha_P_int * num_age_groups;
      for (int j=0; j<alpha_P_int; j++){
        int offset = PStart + j; 
        P[offset] = 0;
      }
      
      int IMStart = i * alpha_IM_int + region * alpha_IM_int * num_age_groups;
      for (int j=0; j<alpha_IM_int; j++){
        int offset = IMStart + j; 
        if (j == 0){
           IM[offset] = compartment_init[1];
          }
         else{
          IM[offset] = 0;
         }
        new_mild += IM[offset];
      }

      int ISStart = i * alpha_IS_int + region * alpha_IS_int * num_age_groups;
      for (int j=0; j<alpha_IS_int; j++){
        int offset = ISStart + j; 
        if (j == 0){
           IS[offset] = compartment_init[2];
          }
         else{
          IS[offset] = 0;
        }
         new_IS += IS[offset];
      }

      int IH1Start = i * alpha_IH1_int + region * alpha_IH1_int * num_age_groups;
      for (int j=0; j<alpha_IH1_int; j++){
        int offset = IH1Start + j; 
        IH1[offset] = 0;
      }

      int IH2Start = i * alpha_IH2_int + region * alpha_IH2_int * num_age_groups;
      for (int j=0; j<alpha_IH2_int; j++){
        int offset = IH2Start + j; 
        IH2[offset] = 0;
      }

      int IH3Start = i * alpha_IH3_int + region * alpha_IH3_int * num_age_groups;
      for (int j=0; j<alpha_IH3_int; j++){
        int offset = IH3Start + j; 
        IH3[offset] = 0;
      }

      int IC2Start = i * alpha_IC2_int + region * alpha_IC2_int * num_age_groups;
      for (int j=0; j<alpha_IC2_int; j++){
        int offset = IC2Start + j; 
        IC2[offset] = 0;
      }

      int IC3Start = i * alpha_IC3_int + region * alpha_IC3_int * num_age_groups;
      for (int j=0; j<alpha_IC3_int; j++){
        int offset = IC3Start + j; 
        IC3[offset] = 0;
      }

      int IH4Start = i * alpha_IH4_int + region * alpha_IH4_int * num_age_groups;
      for (int j=0; j<alpha_IH4_int; j++){
        int offset = IH4Start + j; 
        IH4[offset] = 0;
      }
      
      S[i + region * num_age_groups] = N[i + region * num_age_groups] - new_mild - new_asymp - new_IS;
      R[i + region * num_age_groups] = 0;
      D[i + region * num_age_groups] = 0;
      new_mild_infections[i + region * num_age_groups] = new_mild;
      new_symptomatic_infections[i + region * num_age_groups] = new_IS + new_mild;
      new_deaths[i + region * num_age_groups] = 0;
      new_IH1[i + region * num_age_groups] = 0;
      new_IC2[i + region * num_age_groups] = 0;
      new_IC3[i + region * num_age_groups] = 0;
      new_IH4[i + region * num_age_groups] = 0;
  }

}
