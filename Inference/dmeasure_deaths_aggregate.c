const double *new_deaths = &new_deaths_1_1;
const double *ObsDeaths = &ObsDeaths_1;
const int num_age_groups = n_age_groups;
double lik_total = (give_log) ? 0 : 1;

for (int region=0; region<n_regions; region++){
    double region_lik;

    if (ISNA(ObsDeaths[region])) {
        region_lik = (give_log) ? 0 : 1;
    }
    else{
        double agg_new_D = 0;
        for (int i=0; i<num_age_groups; i++){
            agg_new_D += new_deaths[i + region * num_age_groups];
        }

        int true_ObsDeaths;
        if (t < t_reporting_adjustment){
            double frac_hospitalized_deaths = runif(lower_bound_reporting_uncertainty, 1);
            true_ObsDeaths  = ObsDeaths[region] / (1-(frac_underreported*frac_hospitalized_deaths)); 
        } else{
            true_ObsDeaths = ObsDeaths[region];
        }

        if (true_ObsDeaths <= agg_new_D){
            region_lik = dbetabinom(true_ObsDeaths, agg_new_D, nu_3 * theta_test, dispersion, give_log); 
        } else{
            region_lik = (give_log) ? -1e10 : 0; 
        }

    }

    if (give_log){
        lik_total += region_lik;
    } else{
        lik_total = lik_total * region_lik;
    }
}

lik = lik_total;