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
const double *inv_mu_h = &inv_mu_h_1;

const double *inv_gamma_h = &inv_gamma_h_1;
const double *inv_gamma_h_intercept = &inv_gamma_h_intercept_1;

const double *inv_zeta_icu = &inv_zeta_icu_1;
const double *inv_zeta_icu_intercept = &inv_zeta_icu_intercept_1;
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
    double *IH1 = &IH1_1_1_1;
    double *IH2 = &IH2_1_1_1;

// variant trackers
    double *Evar = &Evar1_1_1;
    double *Pvar = &Pvar1_1_1;
    double *IMvar = &IMvar1_1_1;
    double *IM_deadvar = &IM_deadvar1_1_1;
    double *ISvar = &ISvar1_1_1;

// Tracking states for debugging
    
    // Tracking for states with and without process noise
    double *transmissionRate = &transmissionRate_1_1;
    double *HFRtrack = &HFRtrack_1_1;
    double *ICUtrack = &ICUtrack_1_1;
    double *DeathReportTrack = &DeathReportTrack_1_1;
    double *IHRtrack = &IHRtrack_1_1; 
    double *IFRtrack = &IFRtrack_1_1;
    double *freqtrack = &freqtrack_1_1;

// accumulators
    double *new_deaths = &new_deaths_1_1;
    double *new_hosp_deaths = &new_hosp_deaths_1_1;
    double *new_nonhosp_deaths = &new_nonhosp_deaths_1_1;
    double *new_hospitalizations = &new_hospitalizations_1_1;
    double *Inc = &Inc_1_1;

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



for (int region=0; region<n_regions; region += 1){

    if (t == 384){ // If we're initializing the new variant, move infectious people over to their appropriate compartment
      move_classes(E, Evar, alpha_E_int, init_freq, region);
      move_classes(P, Pvar, alpha_P_int, init_freq, region);
      move_classes(IM, IMvar, alpha_IM_int, init_freq, region);
      move_classes(IM_dead, IM_deadvar, alpha_IM_int, init_freq, region);
      move_classes(IS, ISvar, alpha_IS_int, init_freq, region);
    }

    // Calculate non-hospital death underreporting
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

    double gamma_h = 1 / (inv_gamma_h[region] + inv_gamma_h_intercept[region]);
    double mu_h = 1/ (inv_mu_h[region]);
    double zeta_icu = 1/(inv_zeta_icu[region] + inv_zeta_icu_intercept[region]);
    if (inv_gamma_h[region] + inv_gamma_h_intercept[region] <=0 ){
      gamma_h = 1000;
    }
    if (inv_zeta_icu[region] + inv_zeta_icu_intercept[region] <=0){
      zeta_icu = 1000;
    }
    
    // Figure out total infectious
      double infectious = get_sum(region, alpha_P_int, P) + get_sum(region, alpha_IM_int, IM) + get_sum(region, alpha_IM_int, IM_dead) + get_sum(region, alpha_IS_int, IS);
      double infectious_new = get_sum(region, alpha_P_int, Pvar) + get_sum(region, alpha_IM_int, IMvar) + get_sum(region, alpha_IM_int, IM_deadvar) + get_sum(region, alpha_IS_int, ISvar);
      
      freqtrack[region] = infectious_new / (infectious_new+infectious);

      // Calculate tranmission rate and add in process noise
      double process_noise = sqrt(dt) * eta * rnorm(0,1);
      int n_changepoint_int = n_changepoints[region];
      int changepoint_offset = get_changepoint_index(region, n_changepoints);
      double transmission = get_transmission(changepoint_offset, n_changepoint_int, t);
      transmission += transmission * process_noise;
      transmissionRate[region] = transmission;
      double transmission_var = transmission * scalefactor;

    // Calculate time-varying ratios
    int n_hfr_changepoint_int = n_HFR_changepoints[region];
    int hfr_changepoint_offset = get_changepoint_index(region, n_HFR_changepoints);
    double HFR = calc_time_varying_param(hfr_changepoint_offset, n_hfr_changepoint_int, t, HFR_changepoint_values, HFR_values, HFR_max, HFR_min);
    double HFRvar = HFR * scalefactor_death;
    HFRtrack[region] = HFR;
    int n_icu_changepoint_int = n_ICU_changepoints[region];
    int icu_changepoint_offset = get_changepoint_index(region, n_ICU_changepoints);
    double ICUfrac = calc_time_varying_param(icu_changepoint_offset, n_icu_changepoint_int, t, ICU_changepoint_values, ICU_values, ICU_max, ICU_min);
    double phi = unlogit(phi_scale[region], IFR_constraint, 0);
    double phivar = phi * scalefactor_death;

    double IHR;
    if (t == 47){
      IHR = (IFR_constraint - phi) / (HFR - phi);
    } else{
      IHR = IHRtrack[region];
    }
    IHRtrack[region] = IHR;

    if (IHR < 0){
      IHR = 0;
    }
    IFRtrack[region] = HFR * IHR + (1-IHR) * phi;

      // Outflow from S, consider as competing rates problem
      double lambda = (transmission) * infectious / N[region];
      double lambdavar = (transmission_var) * infectious_new / N[region];
      double dS[2];
      double rates[2] = {lambda, lambdavar};
      reulermultinom(2, S[region], rates, dt, dS);
      S[region] += -dS[0] - dS[1];



 // Outflows from E
      double dE[alpha_E_int];
      double dEvar[alpha_E_int];
      update_states(region, alpha_E_int, dE, E, dt, sigma, dS[0]);
      update_states(region, alpha_E_int, dEvar, Evar, dt, sigma, dS[1]);

      // Outflows from P
      double dP[alpha_P_int];
      double dPvar[alpha_P_int];
      update_states(region, alpha_P_int, dP, P, dt, zeta_s, dE[alpha_E_int-1]);
      update_states(region, alpha_P_int, dPvar, Pvar, dt, zeta_s, dEvar[alpha_E_int-1]);
      
      // Split the last outflow from P into severe and mild
      double dP_s = rbinom(dP[alpha_P_int-1], IHR);
      double dP_m = dP[alpha_P_int-1] - dP_s;
      double dP_m_to_death = rbinom(dP_m, phi); // phi is fraction of non-severe infections that die
      double dP_m_to_recover = dP_m - dP_m_to_death;
      
      // Split the last outflow from P into severe and mild
      double dPvar_s = rbinom(dPvar[alpha_P_int-1], IHR);
      double dPvar_m = dPvar[alpha_P_int-1] - dPvar_s;
      double dPvar_m_to_death = rbinom(dPvar_m, phivar); // phi is fraction of non-severe infections that die
      double dPvar_m_to_recover = dPvar_m - dPvar_m_to_death;

      // Outflows from IM --> Recover
      double dIM[alpha_IM_int];
      update_states(region, alpha_IM_int, dIM, IM, dt, gamma_m, dP_m_to_recover);
      double dIMvar[alpha_IM_int];
      update_states(region, alpha_IM_int, dIMvar, IMvar, dt, gamma_m, dPvar_m_to_recover);

      // Outflows from IM --> Death
      double dIM_dead[alpha_IM_int];
      update_states(region, alpha_IM_int, dIM_dead, IM_dead, dt, mu_m, dP_m_to_death);
      double dIM_deadvar[alpha_IM_int];
      update_states(region, alpha_IM_int, dIM_deadvar, IM_deadvar, dt, mu_m, dPvar_m_to_death);

      // Severe, but not yet hospitalized
      double dIS[alpha_IS_int];
      update_states(region, alpha_IS_int, dIS, IS, dt, zeta_h, dP_s);
      double dISvar[alpha_IS_int];
      update_states(region, alpha_IS_int, dISvar, ISvar, dt, zeta_h, dPvar_s);

      // Split last outflow of IS into: deaths, ICU and recover, and recover
      double S_to_recover = 0;
      double S_to_death = 0;
      double S_to_ICU;
      S_to_recover = rbinom(dIS[alpha_IS_int-1], (1-HFR));
      S_to_death = dIS[alpha_IS_int-1] - S_to_recover;
      S_to_ICU = rbinom(dIS[alpha_IS_int-1], ICUfrac);

      double Svar_to_recover = 0;
      double Svar_to_death = 0;
      double Svar_to_ICU;
      Svar_to_recover = rbinom(dISvar[alpha_IS_int-1], (1-HFRvar));
      Svar_to_death = dISvar[alpha_IS_int-1] - Svar_to_recover;
      Svar_to_ICU = rbinom(dISvar[alpha_IS_int-1], ICUfrac);

      // Path through ICU
      double dpreIC[alpha_preIC_int];
      update_states(region, alpha_preIC_int, dpreIC, preIC, dt, pre_icu_rate, S_to_ICU + Svar_to_ICU);
      double dIC[alpha_IC_int];
      update_states(region, alpha_IC_int, dIC, IC, dt, zeta_icu, dpreIC[alpha_preIC_int - 1]);

      // Hospitalized -> Recover
      double dIH1[alpha_IH1_int];
      update_states(region, alpha_IH1_int, dIH1, IH1, dt, gamma_h, S_to_recover + Svar_to_recover);

      // Hospitalized -> Dead
      double dIH2[alpha_IH2_int];
      update_states(region, alpha_IH2_int, dIH2, IH2, dt, mu_h, S_to_death + Svar_to_death);

      // Update recovered and dead
      R[region] += dIM[alpha_IM_int-1] + dIH1[alpha_IH1_int-1] + dIMvar[alpha_IM_int-1];
      D[region] += dIH2[alpha_IH2_int-1] + dIM_dead[alpha_IM_int-1] + dIM_deadvar[alpha_IM_int-1];
      
      // accumulators
      new_deaths[region] += dIH2[alpha_IH2_int-1] + dIM_dead[alpha_IM_int-1]+ dIM_deadvar[alpha_IM_int-1];
      new_hosp_deaths[region] += dIH2[alpha_IH2_int-1];
      new_nonhosp_deaths[region] += dIM_dead[alpha_IM_int-1]+ dIM_deadvar[alpha_IM_int-1];
      new_hospitalizations[region] += dIS[alpha_IS_int-1] + dISvar[alpha_IS_int-1];
      Inc[region] += dE[alpha_E_int - 1] + dEvar[alpha_E_int - 1];
}
