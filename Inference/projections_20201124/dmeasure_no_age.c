const int reg = round(region_to_test);

const double *new_deaths = &new_deaths_1_1; // latent incident deaths
const double *new_hosp_deaths = &new_hosp_deaths_1_1; // latent hospitalized deaths
const double *new_nonhosp_deaths = &new_nonhosp_deaths_1_1; // latent hospitalized deaths
const double *hosp_death_report = &death_reporting_1;
const double *hosp_report = &hosp_reporting_1;

double *ObsDeaths = &ObsDeaths_1; // observed incident deaths in each region
double *ObsHospDeaths = &ObsHospDeaths_1; // observed incident deaths in each region
double *ObsHosp = &ObsHosp_1; // observed total census hospitalizations in each region
double *IFRtrack = &IFRtrack_1_1; 
const double *IH1 = &IH1_1_1_1;
const double *IH2 = &IH2_1_1_1;

int alpha_IH1_int = round(alpha_IH1);
int alpha_IH2_int = round(alpha_IH2); 

// Start a counter for likelihood
double lik_total = 0;

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

// Set up reporting for deaths, not region-specific
double excess = observed_all_cause - (exp_all_cause + rgeom(pr));
if (excess < 0){
    excess = 0;
}
double total_death_reporting;
if (ll_deaths > excess | excess == 0){
    total_death_reporting = 0.999;
} else{
    total_death_reporting = ll_deaths / excess;
}

for (int region=start_loop; region<end_loop; region += 1){
    if (IFRtrack[region] > 0.01){
        lik_total += -1e3;
    }

    double hosp_death_reporting = hosp_death_report[region];
    double total_hospital_reporting= hosp_report[region];
    if (hosp_death_reporting >= 1){
        hosp_death_reporting = 0.999;
    }
    if (total_hospital_reporting >= 1){
        total_hospital_reporting = 0.999;
    }
    
    // Total deaths
    double region_lik_deaths;
    if (ISNA(ObsDeaths[region])) {
        region_lik_deaths = 0;
    }
    else{
        double agg_new_D = new_deaths[region];
        region_lik_deaths = dpois(ObsDeaths[region], agg_new_D * total_death_reporting, 1); 
    }
    lik_total += region_lik_deaths;

    //Hospital deaths
    double region_lik_hosp_deaths;
    if (ISNA(ObsHospDeaths[region])) {
        region_lik_hosp_deaths = 0;
    }
    else{
        double agg_new_D = new_hosp_deaths[region];
        region_lik_hosp_deaths = dpois(ObsHospDeaths[region], agg_new_D * hosp_death_reporting, 1); 
    }
    lik_total += region_lik_hosp_deaths;

    // Census hospitalization
    double region_lik_hosp;
    if (ISNA(ObsHosp[region])) {
        region_lik_hosp = 0;
    }
    else{
        double agg_new_hosp = 0;
        int IH1_Start = region * alpha_IH1_int;
        for (int k=0; k<alpha_IH1_int; k++){
            int offset = IH1_Start + k;
            agg_new_hosp += IH1[offset];
        }
        int IH2_Start = region * alpha_IH2_int;
        for (int k=0; k<alpha_IH2_int; k++){
            int offset = IH2_Start + k;
            agg_new_hosp += IH2[offset];
        }
        region_lik_hosp = dpois(ObsHosp[region], agg_new_hosp * total_hospital_reporting, 1); 
    }
    lik_total += region_lik_hosp;
}
lik = (give_log) ? lik_total : exp(lik_total);