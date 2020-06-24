const double *new_deaths = &new_deaths_1_1; // latent incident deaths
double *ObsDeaths = &ObsDeaths_1;
double *ObsICU = &ObsICU_1; // observed people in the ICU in each region
const double *IC2 = &IC2_1_1_1; // latent covid-positive ICU
const double *IC3 = &IC3_1_1_1; // latent covid-positive ICU
const int num_age_groups = n_age_groups;
const double *frac_underreported = &frac_underreported_1;
const double *frac_underreported_se = &frac_underreported_se_1;
const double *icu_reporting = &icu_report_1;
int alpha_IC2_int = alpha_IC2; 
int alpha_IC3_int = alpha_IC3;

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
    // loop and aggregate latent incident deaths accross age group
    double agg_new_D = 0;
    
    double frac_underreported_draw = rnorm(frac_underreported[region], frac_underreported_se[region]);
    if (frac_underreported_draw > 1){
        frac_underreported_draw=1;
    } else if(frac_underreported_draw < 0){
        frac_underreported_draw=0;
    }
    double reporting =  1-frac_underreported_draw;
    if (reporting >= 1){
        reporting = runif(0.95, 0.999);
    }

    for (int i=0; i<num_age_groups; i++){
        agg_new_D += new_deaths[i + region * num_age_groups];
    }
    // Draw hospital deaths
    double d = rbinom(agg_new_D, reporting);
    ObsDeaths[region] = d; 
    

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
    double icu_report = rnorm(icu_reporting[region], 0.1);
    if (icu_report > 1){
        icu_report=1;
    } else if(icu_report < 0){
        icu_report=0;
    }
    ObsICU[region] = rbinom(agg_new_ICU, icu_report);
}
