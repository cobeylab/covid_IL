const double *new_hosp_deaths = &new_hosp_deaths_1_1; // latent incident deaths
const double *new_deaths = &new_deaths_1_1; // latent incident deaths
const double *ObsHospDeaths = &ObsHospDeaths_1; // observed incident deaths in each region
const double *ObsNonHospDeaths = &ObsNonHospDeaths_1; // observed incident deaths in each region
const double *ObsICU = &ObsICU_1; // observed people in the ICU in each region
const double *IC2 = &IC2_1_1_1; // latent covid-positive ICU
const double *IC3 = &IC3_1_1_1; // latent covid-positive ICU
const int num_age_groups = n_age_groups;
int alpha_IC2_int = alpha_IC2; 
int alpha_IC3_int = alpha_IC3;
const double *frac_underreported=&frac_underreported_1;
const double *frac_underreported_se = &frac_underreported_se_1;
// Start a counter for likelihood
double lik_total = 0;

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
    double frac_underreported_draw = rnorm(frac_underreported[region], frac_underreported_se[region]);
    if (frac_underreported_draw > 1){
        frac_underreported_draw=1;
    } else if(frac_underreported_draw < 0){
        frac_underreported_draw=0;
    }
    double hosp_reporting_scale;
    if (t < t_reporting_adjustment){
        double frac_hospitalized_deaths = runif(lower_bound_reporting_uncertainty, 0.1);
        hosp_reporting_scale =  (1-(frac_underreported_draw*frac_hospitalized_deaths));
    } else{
        hosp_reporting_scale=1;
    }


    // check if deaths were observed
    double region_lik_hosp_deaths;
    if (ISNA(ObsHospDeaths[region])) {
        region_lik_hosp_deaths = 0;
    }

    else{
        // loop and aggregate latent incident deaths accross age group
        double agg_new_D = 0;
        for (int i=0; i<num_age_groups; i++){
            agg_new_D += new_hosp_deaths[i + region * num_age_groups];
        }

        // Calculate likelihood
        if (ObsHospDeaths[region] <= agg_new_D){
            region_lik_hosp_deaths = dbetabinom(ObsHospDeaths[region], agg_new_D, nu_3 * theta_test * hosp_reporting_scale, dispersion, 1); 
        } else{
            region_lik_hosp_deaths = -1e10; 
        }

    }

    lik_total += region_lik_hosp_deaths;

    //Rprintf("total lik is %f, hosp lik is %f\n", lik_total, region_lik_hosp_deaths);

    //non-hospital deaths
    double region_lik_nonhosp_deaths;
    if (ISNA(ObsNonHospDeaths[region])) {
        region_lik_nonhosp_deaths = 0;
    }

    else{
        // loop and aggregate latent incident deaths accross age group
        double agg_new_D = 0;
        for (int i=0; i<num_age_groups; i++){
            agg_new_D += new_deaths[i + region * num_age_groups] - new_hosp_deaths[i + region * num_age_groups];
        }

        // Calculate likelihood
        if (ObsNonHospDeaths[region] <= agg_new_D){
            double obsprob = nu_m * theta_test * (1-frac_underreported_draw);

            if (obsprob == 0 & ObsNonHospDeaths[region]==0){
                region_lik_nonhosp_deaths = 0;
            } else if (obsprob == 0 & ObsNonHospDeaths[region]>0){
                region_lik_nonhosp_deaths = -1e10;
            } else{
                region_lik_nonhosp_deaths = dbetabinom(ObsNonHospDeaths[region], agg_new_D, obsprob, dispersion, 1); 
            }

        } else{
            region_lik_nonhosp_deaths = -1e10;
        }
                
    }

    lik_total += region_lik_nonhosp_deaths;

    //Rprintf("total lik is %f, nonhosp lik is %f\n", lik_total, region_lik_nonhosp_deaths);
    // check if ICU cases were observed
    double region_lik_ICU;
    if (ISNA(ObsICU[region])) {
        region_lik_ICU = 0;
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
            region_lik_ICU = dbetabinom(ObsICU[region], agg_new_ICU, nu_3 * theta_test * hosp_reporting_scale, dispersion, 1); 
        } else{
            region_lik_ICU = -1e10; 
        }

    }


    lik_total += region_lik_ICU;

    //if (ISNA(lik_total)){
    //    Rprintf("%f %f %f region is %d time is %f\n", region_lik_hosp_deaths, region_lik_nonhosp_deaths, region_lik_ICU, region, t);
    //}

    //if (t==139){
    //    Rprintf("total lik is %f, icu lik is %f, nonhosp lik is %f, hosp lik is %f, region is %d\n", lik_total, region_lik_ICU, region_lik_nonhosp_deaths, region_lik_hosp_deaths, region);
    //}

}


lik = (give_log) ? lik_total : exp(lik_total);