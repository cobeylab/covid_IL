const int num_age_groups = n_age_groups;
const int reg = round(region_to_test);

const double *new_deaths = &new_deaths_1_1; // latent incident deaths
const double *new_hosp_deaths = &new_hosp_deaths_1_1; // latent hospitalized deaths
const double *new_nonhosp_deaths = &new_nonhosp_deaths_1_1; // latent hospitalized deaths
const double *region_non_hosp = &region_non_hosp_1; // region-specific

double *ObsDeaths = &ObsDeaths_1; // observed incident deaths in each region
double *ObsHospDeaths = &ObsHospDeaths_1; // observed incident deaths in each region
double *ObsICU = &ObsICU_1; // observed people in the ICU in each region
double *ObsHosp = &ObsHosp_1; // observed census hospitalizations in each region

const double *IC2 = &IC2_1_1_1; // latent covid-positive ICU
const double *IC3 = &IC3_1_1_1; // latent covid-positive ICU

const double *IH1 = &IH1_1_1_1;
const double *IH2 = &IH2_1_1_1;
const double *IH3 = &IH3_1_1_1;
const double *IH4 = &IH4_1_1_1;

const double *frac_underreported=&frac_underreported_1;
const double *frac_underreported_se = &frac_underreported_se_1;

int alpha_IH1_int = round(alpha_IH1);
int alpha_IH2_int = round(alpha_IH2); 
int alpha_IH3_int = round(alpha_IH3);
int alpha_IH4_int = round(alpha_IH4);
int alpha_IC2_int = round(alpha_IC2); 
int alpha_IC3_int = round(alpha_IC3);


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
    // Set region-specific parameters:

    // loop and aggregate latent incident deaths accross age group
    double agg_new_D = 0;
    double frac_underreported_draw = rnorm(frac_underreported[0], frac_underreported_se[0]);
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
    // Draw total deaths
    double d = rbinom(agg_new_D, reporting);
    ObsDeaths[region] = d; 
    

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
    double icu_reporting = rnorm(0.95, 0.01);
    if (icu_reporting >= 1){
        icu_reporting=0.999;
    } else if(icu_reporting < 0){
        icu_reporting=0;
    }
    ObsICU[region] = rbinom(agg_new_ICU, icu_reporting);

    // aggregate latent hosp deaths
    double agg_new_hospD = 0;
    for (int i=0; i<num_age_groups; i++){
        agg_new_hospD += new_hosp_deaths[i + region * num_age_groups];
    }
    double hosp_death_reporting = rnorm(0.95, 0.01);
    if (hosp_death_reporting >= 1){
        hosp_death_reporting=0.999;
    } else if(hosp_death_reporting < 0){
        hosp_death_reporting=0;
    }
    ObsHospDeaths[region] = rbinom(agg_new_hospD, hosp_death_reporting);      
    

    // Aggregate latent hospitalizations
    double agg_new_hosp = 0;
    for (int i=0; i<num_age_groups; i++){
        int IH1_Start = i * alpha_IH1_int + region * alpha_IH1_int * num_age_groups;
        for (int k=0; k<alpha_IH1_int; k++){
            int offset = IH1_Start + k;
            agg_new_hosp += IH1[offset];
        }
        int IH2_Start = i * alpha_IH2_int + region * alpha_IH2_int * num_age_groups;
        for (int k=0; k<alpha_IH2_int; k++){
            int offset = IH2_Start + k;
            agg_new_hosp += IH2[offset];
        }
        int IH3_Start = i * alpha_IH3_int + region * alpha_IH3_int * num_age_groups;
        for (int k=0; k<alpha_IH3_int; k++){
            int offset = IH3_Start + k;
            agg_new_hosp += IH3[offset];
        }
        int IH4_Start = i * alpha_IH4_int + region * alpha_IH4_int * num_age_groups;
        for (int k=0; k<alpha_IH4_int; k++){
            int offset = IH4_Start + k;
            agg_new_hosp += IH4[offset];
        }
    }
    double hosp_reporting = rnorm(0.95, 0.01);
    if (hosp_reporting >= 1){
        hosp_reporting=0.999;
    } else if(hosp_reporting < 0){
        hosp_reporting=0;
    }
    ObsHosp[region] = rbinom(agg_new_hosp, hosp_reporting);

}
