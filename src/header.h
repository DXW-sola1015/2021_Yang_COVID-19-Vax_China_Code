#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>

#define float double

//GG -- add definitions of global variables

typedef struct
{
  /*from parameters file*/
  double gamma; /* recovery rate from each stage of infection */
    int Nit; /* number of iterations */
  int Tmax; /* number of simulated days */
  int ZETA;
  double DELTAT; /* time step */
  
  /*command line parameters*/
  int is; /*selector of type of susceptibility; 1: age-specific susceptibility; 0: homogeneous susceptibility*/
  int interv; /*identifier of the simulated vaccination scenario*/
  double ni_mean; /*number of daily importations*/
  double omega_1; /*rates of transition between rump up stages*/
  double omega_2; /*rates of transition between rump up stages*/
  double waning_rate; /*1/waning rate duration of protection*/
  double perc_sym_detected;
  int tstartvax; /*time at which vaccination starts*/
  int tstartinf; /*time at which transmission starts*/
  int cov; /*percentage of people to be vaccinated at the end*/
  int capacity; /*number of first doses the system is capable to administer daily. Second doses will be administered with a delay given by the rampup rate*/
  int vax_only_susc; /*identifier of population eligible for vaccination*/
  int vax_prevent; /*identifier of effect of vaccination, 0: prevent infection, 1:prevent symptoms*/
  
  
  double prob_gamma; /*binomial probability of progressing to recovery at a given time step*/
  double prob_omega_1; /*binomial probability of progressing from V0 to V1 at a given time step*/
  double prob_omega_2; /*binomial probability of progressing from V1 to V2 at a given time step*/
  double prob_waning; /*binomial probability of progressing from V2 to W at a given time step*/
  double prob_foi; /* binomial probability of progrssing from S/W to I at a given time step*/
  double prob_foi_0;  /* binomial probability of progressing from V0 to I at a given time step*/
  double prob_foi_1; /* binomial probability of progressing from V1 to I at a given time step*/
  double prob_foi_2; /* binomial probability of progressing from V2 to I at a given time step*/

  
  /*population*/
  int pop; /*total population of the state*/
  int n_ages; /* number of age groups */
  int *age_lim_lo; /*lower limits age groups*/
  int *age_lim_hi; /*upper limits age groups*/
  int *N_u; /*population by age with underlying conditions in the considered state (from file)*/
  int *N_nu; /*population by age without underlying conditions in the considered state (from file)*/
  int *N; /*population by age in the considered state (sum of N_u and N_nu)*/
  /*Other variables*/
  double *beta; /*transmission rate, including transmissibility reduction due to population awareness phi (Eq4)*/
  double *psym; /*probability of developing symptoms (Poletti et al 2020)*/

  double **rel_sus; /*relative susceptibility by age */
  double **prob_immune; /*fraction of immune by age at initialization*/

  /*Contact matrices computed from contact diaries (POLYMOD) (Eq1)*/
  double **CMhouse; /*average contacts at home H*/
  double **CMschool; /*average contacts at school S*/
  double **CMworkall; /*average contacts at work (average over workers and not)*/
  double **CMrandom; /*average contacts in the community T+L+O*/
  /*age-group-specific overall contact matrix*/
  double **CMtot;

  /*vaccination*/
  double **vax_eff; /*vaccine efficacy by age and ramp up stage*/
  int **vax_cov; /*support*/
  int **vax_cov_u; /*number of doses to be administered in each time step to individuals with underlying conditions*/
  int **vax_cov_tmp_u; 
  int **vax_cov_nu; /*number of doses to be administered  in each time step to individuals without underlying conditions*/
  int **vax_cov_tmp_nu; 
} Param;



/*--------------------------------------------------------*/
int get_line(char **line,FILE *fp);
char *get_value(char *cl[],int ncl,char opt[]);

/*Functions to read inputs*/
/*population by age in the state*/
int read_age_structure(char file[],char sep,Param *param);

 /*parameters file*/
int read_param(char file[],char sep, Param *param);

/*transmission rates*/
int read_beta(char file[],char sep,Param *param);

/*initial percentage of immune by age*/
int read_initial_immunity_uncertainty(char file[],char sep,Param *param);

/*age-specific susceptibility to infection*/
int read_relative_susceptibility_uncertainty(char file[],char sep,Param *param); 

 /*contact matrices computed from contact diaries (POLYMOD)*/
int read_contact_matrix_by_setting(char file[],char sep,Param *param,
				   int setting);
/*vaccine efficacy by age and rump up stage*/
int read_vaccine_efficacy(char file[],char sep,Param *param);

/*vaccination coverage by age and over time*/
int read_vaccination_coverage(char file[],char sep,Param *param);

/*combine contact matrices in the different settings to reconstruct matrix of overall contacts*/
int build_total_matrix(Param *param);
