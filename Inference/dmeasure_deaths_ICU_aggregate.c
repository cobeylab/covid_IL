const double *new_hosp_deaths = &new_hosp_deaths_1_1; // latent incident deaths
const double *ObsDeaths = &ObsDeaths_1; // observed incident deaths in each region
const double *ObsICU = &ObsICU_1; // observed people in the ICU in each region
const double *IC2 = &IC2_1_1_1; // latent covid-positive ICU
const double *IC3 = &IC3_1_1_1; // latent covid-positive ICU
const int num_age_groups = n_age_groups;
int alpha_IC2_int = alpha_IC2; 
int alpha_IC3_int = alpha_IC3;
// Start a counter for likelihood
double lik_total = (give_log) ? 0 : 1;

// loop over each region
for (int region=0; region<n_regions; region++){
    // check if deaths were observed
    double region_lik_deaths;
    if (ISNA(ObsDeaths[region])) {
        region_lik_deaths = (give_log) ? 0 : 1;
    }

    else{
        // loop and aggregate latent incident deaths accross age group
        double agg_new_D = 0;
        for (int i=0; i<num_age_groups; i++){
            agg_new_D += new_hosp_deaths[i + region * num_age_groups];
        }
        double frac_hospitalized_deaths = runif(lower_bound_reporting_uncertainty, 1);
        double underreporting =  (1-(frac_underreported*frac_hospitalized_deaths)); 

        // Calculate likelihood
        if (ObsDeaths[region] <= agg_new_D){
            region_lik_deaths = dbetabinom(ObsDeaths[region], agg_new_D, nu_3 * theta_test * underreporting, dispersion, give_log); 
        } else{
            region_lik_deaths = (give_log) ? -1e10 : 0; 
        }

    }

    if (give_log){
        lik_total += region_lik_deaths;
    } else{
        lik_total = lik_total * region_lik_deaths;
    }


    // check if ICU cases were observed
    double region_lik_ICU;
    if (ISNA(ObsICU[region])) {
        region_lik_ICU = (give_log) ? 0 : 1;
    }
    else{

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
        if (ObsICU[region] <= agg_new_ICU){
            region_lik_ICU = dbetabinom(ObsICU[region], agg_new_ICU, nu_3 * theta_test, dispersion, give_log); 
        } else{
            region_lik_ICU = (give_log) ? -1e10 : 0; 
        }

    }

    if (give_log){
        lik_total += region_lik_ICU;
    } else{
        lik_total = lik_total * region_lik_ICU;
    }
}

lik = lik_total;