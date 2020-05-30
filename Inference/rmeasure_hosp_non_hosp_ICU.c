const double *new_hosp_deaths = &new_hosp_deaths_1_1; // latent incident deaths
const double *new_deaths = &new_deaths_1_1; // latent incident deaths
double *ObsHospDeaths = &ObsHospDeaths_1; // observed incident deaths in each region
double *ObsNonHospDeaths = &ObsNonHospDeaths_1; // observed incident deaths in each region
double *ObsICU = &ObsICU_1; // observed people in the ICU in each region
const double *IC2 = &IC2_1_1_1; // latent covid-positive ICU
const double *IC3 = &IC3_1_1_1; // latent covid-positive ICU
const int num_age_groups = n_age_groups;
const double *frac_underreported = &frac_underreported_1;
int alpha_IC2_int = alpha_IC2; 
int alpha_IC3_int = alpha_IC3;



// loop over each region
for (int region=0; region<n_regions; region++){
    // loop and aggregate latent incident deaths accross age group
    double agg_new_D = 0;
    for (int i=0; i<num_age_groups; i++){
        agg_new_D += new_hosp_deaths[i + region * num_age_groups];
    }

    double reporting;
    if (t < t_reporting_adjustment){
        double frac_hospitalized_deaths = runif(lower_bound_reporting_uncertainty, 0.5);
        reporting =  (1-(frac_underreported[region]*frac_hospitalized_deaths));
    } else{
        reporting=1;
    }
    // Draw hospital deaths
    ObsHospDeaths[region] = rbetabinom(agg_new_D, nu_3 * theta_test * reporting, dispersion); 
    

    //non-hospital deaths
    // loop and aggregate latent incident deaths accross age group
    agg_new_D = 0;
    for (int i=0; i<num_age_groups; i++){
        agg_new_D += new_deaths[i + region * num_age_groups] - new_hosp_deaths[i + region * num_age_groups];
    }
    reporting = 1-frac_underreported[region];
    // Draw nonhosp deaths
    ObsNonHospDeaths[region] = rbetabinom(agg_new_D, nu_m * theta_test * reporting, dispersion);
    // check if ICU cases were observed

    // aggregate latent ICU over all ages and subcompartments
    double agg_new_ICU = 0;
    for (int i=0; i<num_age_groups; i++){
        int IC2_Start = i * alpha_IC2_int + region * alpha_IC2_int * num_age_groups;
        for (int k=0; k<alpha_IC2_int; k++){
            int offset = IC2_Start + k;
            agg_new_ICU += IC2[offset];
        }
        int IC3_Start = i * alpha_IC3_int + region * alpha_IC3_int * num_age_groups;
        for (int k=0; k<alpha_IC3_int; k++){
            int offset = IC3_Start + k;
            agg_new_ICU += IC3[offset];
        }
    }
    // Calculate likelihood
    ObsICU[region] = rbetabinom(agg_new_ICU, nu_3 * theta_test, dispersion);
}
