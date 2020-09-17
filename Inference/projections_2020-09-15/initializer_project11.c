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

// Accumulators
double *new_deaths = &new_deaths_1_1;
double *new_hosp_deaths = &new_hosp_deaths_1_1;
double *new_nonhosp_deaths = &new_nonhosp_deaths_1_1;
double *new_hospitalizations = &new_hospitalizations_1_1;
double *new_symptomatic_infections = &new_symptomatic_infections_1_1;
double *Inc = &Inc_1_1;
//double *new_mild_infections = &new_mild_infections_1_1;
//double *new_IH1 = &new_IH1_1_1;
//double *new_IC2 = &new_IC2_1_1;
//double *new_IC3 = &new_IC3_1_1;
//double *new_IH4 = &new_IH4_1_1;

// initialize parameter pointers
const double *N = &N_1_1;
const double *age_dist = &age_dist_1_1;
const double *num_init = &num_init_1;


//double *nonicubed = &nonicubed_1;
//double *icubed = &icubed_1;

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
int regparam = 0;
for (int region=start_loop; region<end_loop; region += 1 ){

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

  // initialize multinomial draw to determine age distribution
  int rN[num_age_groups];
  int total_pop = 0;
  double age_dist_region[num_age_groups];
  for (int ag=0; ag<num_age_groups; ag+=1){
    age_dist_region[ag] = age_dist[ag + regparam * num_age_groups];
    total_pop += N[ag + region * num_age_groups];
  }
  double num_init_final = num_init[region] * total_pop;
  rmultinom(ceil(num_init_final), age_dist_region, num_age_groups, rN);

  for (int i=0; i<num_age_groups; i += 1){
      int EStart = i * alpha_E_int + region * alpha_E_int * num_age_groups;
      for (int j=0; j<alpha_E_int; j++){
        int offset = EStart + j;
        if (j == 0){
          E[offset] = rN[i];
        } 
        else{
          E[offset] = 0;    
        }

      }

      int AStart = i * alpha_A_int + region * alpha_A_int * num_age_groups;
      for (int j=0; j<alpha_A_int; j++){
        int offset = AStart + j; 
        A[offset] = 0;
      }

      int PStart = i * alpha_P_int + region * alpha_P_int * num_age_groups;
      for (int j=0; j<alpha_P_int; j++){
        int offset = PStart + j; 
        P[offset] = 0;
      }
      
      // divide mild infections into recovery and death pathways
      int IMStart = i * alpha_IM_int + region * alpha_IM_int * num_age_groups;
      for (int j=0; j<alpha_IM_int; j++){
        int offset = IMStart + j; 
        IM[offset] = 0;
        IM_dead[offset] = 0;
      }

      int ISStart = i * alpha_IS_int + region * alpha_IS_int * num_age_groups;
      for (int j=0; j<alpha_IS_int; j++){
        int offset = ISStart + j; 
        IS[offset] = 0;
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
      S[i + region * num_age_groups] = N[i + region * num_age_groups] - rN[i];
      R[i + region * num_age_groups] = 0;
      D[i + region * num_age_groups] = 0;
      new_deaths[i + region * num_age_groups] = 0;
      new_hosp_deaths[i + region * num_age_groups] = 0;
      new_nonhosp_deaths[i + region * num_age_groups] = 0;
      new_hospitalizations[i + region * num_age_groups] = 0;
      Inc[i + region * num_age_groups] = rN[i];
      //new_IH1[i + region * num_age_groups] = 0;
      //new_IC2[i + region * num_age_groups] = 0;
      //new_IC3[i + region * num_age_groups] = 0;
      //new_IH4[i + region * num_age_groups] = 0;
      new_symptomatic_infections[i + region * num_age_groups] = 0;
      //new_mild_infections[i + region * num_age_groups] = 0;

      //if (region == 2){
      //  icubed[i] = 0;
      //  nonicubed[i] = 0;
     //}

  }
}
