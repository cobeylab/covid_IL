// Number of age groups 
  int reg = round(region_to_test);

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

// region-specific fraction of non-hospitalized deaths
  const double *phi_scale = &phi_scale_1;
  const double *inv_mu_0 = &inv_mu_0_1;
  const double *inv_mu_f = &inv_mu_f_1;
  const double *inv_gamma_0 = &inv_gamma_0_1;
  const double *inv_gamma_f = &inv_gamma_f_1;
// region-specific transmission baseline parameter for scaling
  const double *changepoint_values = &changepoints_1_1;
  const double *n_changepoints = &n_changepoints_1;
  const double *beta_values = &beta_values_1_1;

// region-specific HFR
  const double *HFR_changepoint_values = &HFR_changepoint_values_1_1;
  const double *n_HFR_changepoints = &n_HFR_changepoints_1;
  const double *HFR_values = &HFR_values_1_1;

// region-specific fraction ICU
  const double *ICU_changepoint_values = &ICU_changepoint_values_1_1;
  const double *n_ICU_changepoints = &n_ICU_changepoints_1;
  const double *ICU_values = &ICU_values_1_1;
// ICU rate
  const double *inv_zeta_icu_0 = &inv_zeta_icu_0_1;
  const double *inv_zeta_icu_f = &inv_zeta_icu_f_1;

// initialize compartment pointers
    double *S = &S_1_1;
    double *R = &R_1_1;
    double *D = &D_1_1;
    double *E = &E1_1_1;
    double *P = &P1_1_1;
    double *IM = &IM1_1_1;
    double *IM_dead = &IM_dead1_1_1;
    double *IS = &IS1_1_1;
    double *preIC = &preIC_1_1_1;
    double *IC = &IC_1_1_1;

// Chains of severe infectious states
    double *IH1 = &IH1_1_1_1;
    double *IH2 = &IH2_1_1_1;
// Tracking states for debugging
    double *transmissionRate = &transmissionRate_1_1; 
    double *HFRtrack = &HFRtrack_1_1;
    double *phitrack = &phitrack_1_1;
    double *IHRtrack = &IHRtrack_1_1; 
    double *IFRtrack = &IFRtrack_1_1; 
// accumulators
    double *new_deaths = &new_deaths_1_1;
    double *new_hosp_deaths = &new_hosp_deaths_1_1;
    double *new_nonhosp_deaths = &new_nonhosp_deaths_1_1;
    double *new_hospitalizations = &new_hospitalizations_1_1;
    double *Inc = &Inc_1_1;

//Function to make a line
double get_par_value( double p0, double pf, double t0, double tf, double t_now ) {
   double slope = (pf - p0) / (tf - t0);
   double value;
   if ( t_now >= tf){
       value = pf;
   } else if (t_now <= t0){
       value = p0;
   } else{
       value = p0 + slope * (t_now - t0);
   }
  return value;
}
// Function untransform a logit parameter
double unlogit( double logit_value, double pmax, double pmin) {
  return (logit_value * (pmax - pmin)) + pmin;
}
// Function to calculate transmission given a list of changepoints and beta values at those changepoints
double get_transmission( int coffset, int n_timepoints, double t_now){
  double final_beta;
  // if we're at a point before the first timepoint, just assign to first beta value
  if (t_now <= changepoint_values[coffset]){
    final_beta = unlogit(beta_values[coffset], 0.9, 0.4); //constraint on initial beta
  } else if(t_now >= changepoint_values[coffset + n_timepoints - 1]){
    final_beta = beta_values[coffset + n_timepoints - 1];
  } else{
    for (int i=coffset + 1; i <= coffset + n_timepoints - 1; i++){
      if(t_now < changepoint_values[i]){
        double b1 = beta_values[i-1];

        if (i == coffset + 1){
            b1 = unlogit(beta_values[coffset], 0.9, 0.4); //constraint on initial beta
        }

        double b2 = beta_values[i];
        double t1 = changepoint_values[i-1];
        double t2 = changepoint_values[i];
        final_beta = get_par_value(b1, b2, t1, t2, t_now);
        break;
      }
    }
  }
  return final_beta;
}
// Function to calculate transmission given a list of changepoints and HFR values at those changepoints
double get_HFR( int coffset, int n_timepoints, double t_now){
  double final_HFR;
  // if we're at a point before the first timepoint, just assign to first HFR value
  if (t_now <= HFR_changepoint_values[coffset]){
    final_HFR = unlogit(HFR_values[coffset], HFR_max, HFR_min);
  } else if(t_now >= HFR_changepoint_values[coffset + n_timepoints - 1]){
    final_HFR = unlogit(HFR_values[coffset + n_timepoints - 1], HFR_max, HFR_min);
  } else{
    for (int i=coffset + 1; i <= coffset + n_timepoints - 1; i++){
      if(t_now < HFR_changepoint_values[i]){
        double b1 = unlogit(HFR_values[i-1], HFR_max, HFR_min);
        double b2 = unlogit(HFR_values[i], HFR_max, HFR_min);
        double t1 = HFR_changepoint_values[i-1];
        double t2 = HFR_changepoint_values[i];
        final_HFR = get_par_value(b1, b2, t1, t2, t_now);
        break;
      }
    }
  }
  return final_HFR;
}
double get_ICU( int coffset, int n_timepoints, double t_now){
  double final_ICU;
  // if we're at a point before the first timepoint, just assign to first ICU value
  if (t_now <= ICU_changepoint_values[coffset]){
    final_ICU = unlogit(ICU_values[coffset], ICU_max, ICU_min);
  } else if(t_now >= ICU_changepoint_values[coffset + n_timepoints - 1]){
    final_ICU = unlogit(ICU_values[coffset + n_timepoints - 1], ICU_max, ICU_min);
  } else{
    for (int i=coffset + 1; i <= coffset + n_timepoints - 1; i++){
      if(t_now < ICU_changepoint_values[i]){
        double b1 = unlogit(ICU_values[i-1], ICU_max, ICU_min);
        double b2 = unlogit(ICU_values[i], ICU_max, ICU_min);
        double t1 = ICU_changepoint_values[i-1];
        double t2 = ICU_changepoint_values[i];
        final_ICU = get_par_value(b1, b2, t1, t2, t_now);
        break;
      }
    }
  }
  return final_ICU;
}
int get_changepoint_index(int reg){
  int start = 0;
  for (int r=1; r<=reg; r++){
      int z = round(n_changepoints[r-1]);
      start = round(start + z);
  }
  return start;
}
int get_hfr_changepoint_index(int reg){
  int start = 0;
  for (int r=1; r<=reg; r++){
      int z = round(n_HFR_changepoints[r-1]);
      start = round(start + z);
  }
  return start;
}
int get_icu_changepoint_index(int reg){
  int start = 0;
  for (int r=1; r<=reg; r++){
      int z = round(n_ICU_changepoints[r-1]);
      start = round(start + z);
  }
  return start;
}

double HFR;

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

for (int region=start_loop; region<end_loop; region += 1){

    double gamma_h = 1/get_par_value(unlogit(inv_gamma_0[region], gamma_max, gamma_min), unlogit(inv_gamma_f[region], gamma_max, gamma_min), hosp_t_min, hosp_t_max, t); 
    double mu_h = 1/get_par_value(unlogit(inv_mu_0[region], mu_max, mu_min), unlogit(inv_mu_f[region], mu_max, mu_min), hosp_t_min, hosp_t_max, t); 
    double zeta_icu = 1/get_par_value(unlogit(inv_zeta_icu_0[region], zeta_icu_max, zeta_icu_min), unlogit(inv_zeta_icu_f[region], zeta_icu_max, zeta_icu_min), hosp_t_min, hosp_t_max, t); 
    
    double lambda;

      if (region_specific_HFR == 1){
        int n_hfr_changepoint_int = n_HFR_changepoints[region];
        int hfr_changepoint_offset = get_hfr_changepoint_index(region);
        HFR = get_HFR(hfr_changepoint_offset, n_hfr_changepoint_int, t);
      }
      double phi = unlogit(phi_scale[region], IFR_constraint, 0);
      double IHR;
      if (t == 47){
        IHR = (IFR_constraint - phi) / (HFR - phi);
      } else{
        IHR = IHRtrack[region];
      }
      
      HFRtrack[region] = HFR;
      IHRtrack[region] = IHR;

      if (IHR < 0){
        IHR = 0;
      }

    // Figure out fraction going to the ICU
      int n_icu_changepoint_int = n_ICU_changepoints[region];
      int icu_changepoint_offset = get_icu_changepoint_index(region);
      double ICUfrac = get_ICU(icu_changepoint_offset, n_icu_changepoint_int, t);

    // Figure out total infectious
      double infectious = 0;
      int PStart = region * alpha_P_int;
      for (int x=0; x<alpha_P_int; x++){
        int offset = PStart + x;
        infectious += P[offset];
      }
      int IMStart = region * alpha_IM_int;
      for (int x=0; x<alpha_IM_int; x++){
        int offset = IMStart + x;
        infectious += IM[offset] + IM_dead[offset];
      }
      int ISStart = region * alpha_IS_int;
      for (int x=0; x<alpha_IS_int; x++){
        int offset = ISStart + x;
        infectious += IS[offset];
      }

      // Calculate FOI
      double transmission;
      double process_noise = sqrt(dt) * eta * rnorm(0,1);
      int n_changepoint_int = n_changepoints[region];
      int changepoint_offset = get_changepoint_index(region);
      transmission = get_transmission(changepoint_offset, n_changepoint_int, t);
      lambda = (transmission) * infectious / N[region];
      // Add in process noise
      lambda += lambda * process_noise;
      transmissionRate[region] = transmission; //Update the transmission rate state to keep track of transmission

      // Outflow from S
      double dS = rbinom(S[region], 1-exp(-lambda * dt));
      S[region] += -dS;

      // Outflows from E
      int EStart = region * alpha_E_int;
      double dE[alpha_E_int];
      double Pr_E = 1-exp(-sigma * alpha_E * dt);
      dE[0] = rbinom(E[EStart], Pr_E);
      E[EStart] += dS - dE[0];
      
      for (int i=1; i<alpha_E_int; i++){
        int offset = EStart + i;
        dE[i] = rbinom(E[offset], Pr_E);
        E[offset] += dE[i - 1] - dE[i];         
      }
      
      // Outflows from P
      double dP[alpha_P_int];
      double Pr_P = 1-exp(-zeta_s * alpha_P * dt);
      dP[0] = rbinom(P[PStart], Pr_P);
      P[PStart] += dE[alpha_E_int-1] - dP[0];
      for (int i=1; i<alpha_P_int; i++){
        int offset = PStart + i;
        dP[i] = rbinom(P[offset], Pr_P);
        P[offset] += dP[i-1] - dP[i];
      }
      // Split the last outflow from P into severe and mild
      double dP_s = rbinom(dP[alpha_P_int-1], IHR);
      double dP_m = dP[alpha_P_int-1] - dP_s;

      //Rprintf("%f, IFR is %f\n", t, HFR * IHR + (1-IHR) * phi);

      phitrack[region] = (phi * (1-IHR)) / (phi * (1-IHR) + IHR * HFR);
      IFRtrack[region] = HFR * IHR + (1-IHR) * phi;
      
      double dP_m_to_death = rbinom(dP_m, phi); // phi is fraction of non-severe infections that die
      double dP_m_to_recover = rbinom(dP_m, (1-phi));
      
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
      int IM_dead_Start = region * alpha_IM_int;
      double dIM_dead[alpha_IM_int];
      double Pr_IM_dead = 1-exp(-mu_m * alpha_IM * dt); // mu_m is death rate for non-hospitalized deaths
      dIM_dead[0] = rbinom(IM_dead[IM_dead_Start], Pr_IM_dead);
      IM_dead[IM_dead_Start] += dP_m_to_death - dIM_dead[0];
      for (int i=1; i<alpha_IM_int; i++){
        int offset = IM_dead_Start + i;
        dIM_dead[i] = rbinom(IM_dead[offset], Pr_IM_dead);
        IM_dead[offset] += dIM_dead[i - 1] - dIM_dead[i];
      }
      
      // Severe, but not yet hospitalized
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
      double S_to_recover = 0;
      double S_to_death = 0;
      double S_to_ICU;
      S_to_recover = rbinom(dIS[alpha_IS_int-1], (1-HFR));
      S_to_death = dIS[alpha_IS_int-1] - S_to_recover;
      S_to_ICU = rbinom(dIS[alpha_IS_int-1], ICUfrac);


      int preIC_Start = region * alpha_preIC_int;
      double dpreIC[alpha_preIC_int];
      double Pr_preIC = 1-exp(-pre_icu_rate * alpha_preIC_int * dt);
      dpreIC[0] = rbinom(preIC[preIC_Start], Pr_preIC);
      preIC[preIC_Start] += S_to_ICU - dpreIC[0];
      for (int i=1; i<alpha_preIC_int; i++){
        int offset = preIC_Start + i;
        dpreIC[i] = rbinom(preIC[offset], Pr_preIC);
        preIC[offset] += dpreIC[i - 1] - dpreIC[i];
      }

      int IC_Start = region * alpha_IC_int;
      double dIC[alpha_IC_int];
      double Pr_IC = 1-exp(-zeta_icu * alpha_IC_int * dt);

      dIC[0] = rbinom(IC[IC_Start], Pr_IC);

      IC[IC_Start] += dpreIC[alpha_preIC_int - 1] - dIC[0];
      for (int i=1; i<alpha_IC_int; i++){
        int offset = IC_Start + i;
        dIC[i] = rbinom(IC[offset], Pr_IC);
        IC[offset] += dIC[i - 1] - dIC[i];
      }
      


      // Hospitalized -> Recover
      // Outflows from IH1
      int IH1_Start = region * alpha_IH1_int;
      double dIH1[alpha_IH1_int];
      double Pr_IH1 = 1-exp(-gamma_h * alpha_IH1 * dt);
      dIH1[0] = rbinom(IH1[IH1_Start], Pr_IH1);
      IH1[IH1_Start] += S_to_recover - dIH1[0];
      for (int i=1; i<alpha_IH1_int; i++){
        int offset = IH1_Start + i;
        dIH1[i] = rbinom(IH1[offset], Pr_IH1);
        IH1[offset] += dIH1[i - 1] - dIH1[i];
      }
      
      // Hospitalized -> Dead
      // Outflows from IH2
      int IH2_Start = region * alpha_IH2_int;
      double dIH2[alpha_IH2_int];
      double Pr_IH2 = 1-exp(-mu_h * alpha_IH2 * dt);
      dIH2[0] = rbinom(IH2[IH2_Start], Pr_IH2);
      IH2[IH2_Start] += S_to_death - dIH2[0];
      for (int i=1; i<alpha_IH2_int; i++){
        int offset = IH2_Start + i;
        dIH2[i] = rbinom(IH2[offset], Pr_IH2);
        IH2[offset] += dIH2[i - 1] - dIH2[i];
      }

      // Update recovered and dead
      R[region] += dIM[alpha_IM_int-1] + dIH1[alpha_IH1_int-1];
      D[region] += dIH2[alpha_IH2_int-1] + dIM_dead[alpha_IM_int-1];
      
      // accumulators
      new_deaths[region] += dIH2[alpha_IH2_int-1] + dIM_dead[alpha_IM_int-1];
      new_hosp_deaths[region] += dIH2[alpha_IH2_int-1];
      new_nonhosp_deaths[region] += dIM_dead[alpha_IM_int-1];
      new_hospitalizations[region] += dIS[alpha_IS_int-1];
      Inc[region] += dE[alpha_E_int - 1];
}
