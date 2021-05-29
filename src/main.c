#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "header.h"
#include "global.h"
#include <math.h>
#include <sys/stat.h>

int main(int argc, char *argv[]){

  char help[]="gen\n exp_dir\n -w [country]\n -s [seed]\n -selReff [selector of transmissibility]\n -selVAXEFF [selector for vaccine efficacy file, indicates vaccine efficacy achieved after the 2nd dose in <60y]\n -selIMM [identifier of % of immune at initialization]\n -is [select susceptibility by age: 0=homogeneous susceptibility, 1=age-specific susceptibility]\n -cov [vaccination coverage - used to select input file of number of doses administered over time]\n -capacity [system capacity - daily number of first doses the system is capable to administer - used to select input file of number of doses administered over time]\n -interv [identifier of population priority]\n -ni [initial number of infections]\n  -omega_1 [rate of transition between V0 and V1 (interval between 1st and 2nd dose)]\n -omega_2 [rate of transition between V1 and V2 (time needed by 2nd dose to become effective)]\n -waning_rate [rate of transition between V2 and W] - tstartvax [time at which vaccination starts] - tstartinf [time at which trabsmission starts]\n -pdetected [percentage of symptomatic cases detected]\n -vax_only_susc [0=vaccination to susceptibles and undetected 1=vaccination only to susceptibles, i.e. set psym=1 for all ages and pdetected=1] \n -vax_prevent [0=infection 1=symptoms] \n"; 

  /*selectors*/
  char *char_seed=NULL;
  char *char_ni=NULL;
  char *char_where=NULL;
  char *char_phi=NULL;
  char *char_perc_sym_detected=NULL;
  char *char_vaxeff=NULL;
  char *char_iscen=NULL; 
  char *char_is_option=NULL;
  char *char_interv_option=NULL;
  char *char_tstartvax_option=NULL;
  char *char_tstartinf_option=NULL;
  char *char_capacity_option=NULL;
  char *char_cov_option=NULL;
  char *char_vax_only_susc=NULL;
  char *char_vax_prevent=NULL;

   
  char char_rs[40];
  char param_file[1000];
  char age_file[1000];
  char cm_file[1000];
  char rs_file[1000];
  char outdir[1000];
  char output_file[1000];
  
  char i_file[1000];
  FILE *input_file;
 
  register int i,j,q;

  int *S_u, *S_nu; /*Susceptible with/without underlying conditions (u.c.)*/
  int *I_u, *I_nu ; /*Infectious with/without u.c.*/
  int *V0_u, *V0_nu ; /*Vaccinated with the 1st dose w/wo u.c.*/
  int *V1_u, *V1_nu; /*Vaccinated with the 2nd dose w/wo u.c. (2nd dose not effective)*/
  int *V2_u, *V2_nu; /*Vaccinated with the 2nd dose w/wo u.c. (2nd dose effective)*/
  int *W_u, *W_nu; /* Vaccinated with two doses w/wo u.c. for whom vaccine-induced immunity has waned*/
  int *R_u, *R_nu; /* Recovered from infection  w/wo u.c.*/
  
  /*New transitions*/
  int *nI_u, *nI_nu; /* new transitions between S and I */
  int *nIdet_u, *nIdet_nu; /*new transitions between S and I detected by the system*/
  int *nI_0_u, *nI_0_nu;  /*new transitions between V0 and I */
  int *nI_1_u, *nI_1_nu; /*new transitions between V1 and I */
  int *nI_2_u, *nI_2_nu; /*new transitions between V2 and I */
  int *nI_3_u, *nI_3_nu; /*new transitions between W and I */ 
  
  int *nV0_u, *nV0_nu; /* new transitions between S and V0 */
  int *nV1_u, *nV1_nu; /* new transitions between V0 and V1 */
  int *nV2_u, *nV2_nu; /* new transitions between V1 and V2 */
  int *nW_u, *nW_nu; /* new transitions between V2 and W */
  int *nR_u, *nR_nu; /* new transitions between I and R */

  /*Counters*/
  int *C_u,*C_nu; /*counter daily number of new infections among people w/wo u.c.*/
  int *C_u_vax,*C_nu_vax; /*counter daily number of new infections among people w/wo u.c. effectively vaccinated (compartment V2)*/
  int *nsecond_doses_u;/*daily counter second doses*/
  int *nsecond_doses_nu;/*daily counter second doses*/
  int **second_doses_over_time;/*counter second doses (daily)*/
  int *wasted_capacity; /*tmp variable*/
  int *wasted_capacity_u; /*tmp variable*/
  int *wasted_capacity_nu; /*tmp variable*/
  
  int nvax_u, nvax_nu; /*counter new first doses administered to susceptibles */
  int  nvax_wasted_u, nvax_wasted_nu;/*counter new first doses administered to undetected infected/removed individuals */

  /*used for checks*/
  int *cumdet_u; /*counter detected infections by age with u.c. */
  int *cumdet_nu; /*counter detected infections by age without u.c.*/
  int *cuminf; /*counter infections by age among unvaccinated*/

  
  int *cumundet_u; /* cumulative undetected infections among susceptibles with u.c.*/
  int *cumundet_nu; /* cumulative undetected infections among susceptibles without u.c.*/


  
  double *foi; /*age-specific force of infection (lambda)*/
  int initial_infections;
  int Itot=0;
  Param param;

  int t;
  int T;

  FILE *fpout;
  FILE *fpout_vax;

 
  /*-----------------------------------------------------*/
  /*             read command line parameters            */
  /*-----------------------------------------------------*/
  /*check number of arguments and  eventually print help*/
  if(argc<38){ 
    fprintf(stderr,"%s\n",help);
    exit(1);
  }


  if((strcmp(argv[1],"-h")==0 || strcmp(argv[1],"--help")==0 || 
		 strcmp(argv[1],"-help")==0)){
    fprintf(stderr,"%s\n",help);
    exit(1);
  }

  /*Experiment directory (first argument command line)*/
  EXP_DIR=argv[1];

  /* Set seed random generator */ 
  if((char_seed=get_value(argv,argc,"-s"))==NULL){
    fprintf(stderr,"Specify seed!\n");
    exit(1);
  }

  fprintf(stderr,"SEED=%s\n",char_seed);
  
  /* number of initial infections */ 
  if((char_ni=get_value(argv,argc,"-ni"))==NULL){
    fprintf(stderr,"Specify number of initial infections!\n");
    exit(1);
  }else{
    sscanf(char_ni,"%lf",&param.ni_mean);
  }
  initial_infections=(int) round(param.ni_mean);
  fprintf(stderr,"INITIAL INFECTIONS=%d\n",initial_infections); 
    

  /* selector STATE */
  if((char_where=get_value(argv,argc,"-w"))==NULL){
    fprintf(stderr,"Specify where!\n");
    exit(1);
  }
  fprintf(stderr, "%s ", char_where);

  /* selector vaccine efficacy file to simulate */
  if((char_vaxeff=get_value(argv,argc,"-selVAXEFF"))==NULL){
    fprintf(stderr,"Specify vaxeff selector!\n");
    exit(1);
  }
  fprintf(stderr, "Efficacy after 2 doses=%s\n", char_vaxeff);

  /*Selector of effective reproduction number to simulate*/ 
  if((char_phi=get_value(argv,argc,"-selReff"))==NULL){
    fprintf(stderr,"Specify selector Reff!\n");
    exit(1);
  }
  fprintf(stderr, "Reff=%s\n", char_phi);

  /*pdetected - percentage of symptomatic cases detected*/ 
  if((char_perc_sym_detected=get_value(argv,argc,"-pdetected"))==NULL)
    param.perc_sym_detected=1.; /* default */
  else
    sscanf(char_perc_sym_detected,"%lf",&param.perc_sym_detected);
  fprintf(stderr,"pdetected=%lf\n",param.perc_sym_detected);

  /*vax_only_susc - identifier population eligible for vaccination*/ 
  if((char_vax_only_susc=get_value(argv,argc,"-vax_only_susc"))==NULL){
    fprintf(stderr,"Specify parameter vax_only_susc!\n");
    exit(1);
  }else{
    sscanf(char_vax_only_susc,"%d",&param.vax_only_susc);
  }
  fprintf(stderr,"vax_only_susc=%d\n",param.vax_only_susc);

   /*vax_prevent - identifier of effect of vaccination: 0 prevent only symptoms, 1 prevent also infection*/ 
  if((char_vax_prevent=get_value(argv,argc,"-vax_prevent"))==NULL){
     fprintf(stderr,"Specify parameter vax_prevent!\n");
     exit(1);
  }else{
    sscanf(char_vax_prevent,"%d",&param.vax_prevent);
  }
  fprintf(stderr,"vax_prevent=%d\n",param.vax_prevent);

  /*identifier of the assumed initial immunity scenario*/
  if((char_iscen=get_value(argv,argc,"-selIMM"))==NULL){
    fprintf(stderr,"Specify immunity!\n");
    exit(1);
  }
  fprintf(stderr, "Immunity=%s\n", char_iscen);

  /*Select type of age-specific susceptibility*/
  if((char_is_option=get_value(argv,argc,"-is"))==NULL){
    fprintf(stderr,"Specify parameter is!\n");
    exit(1);
  }else{
    sscanf(char_is_option,"%d",&param.is);
  }
  fprintf(stderr,"susc=%d\n",param.is);
 
  /*identifier of the priority population*/
  if((char_interv_option=get_value(argv,argc,"-interv"))==NULL){
    fprintf(stderr,"Specify identifier of prioirity population!\n");
    exit(1);
  }else{
    sscanf(char_interv_option,"%d",&param.interv);
  }
  fprintf(stderr,"PRIORITY=%d\n",param.interv);

  /*time at which vaccination starts*/
  if((char_tstartvax_option=get_value(argv,argc,"-tstartvax"))==NULL){
    fprintf(stderr,"Specify day of vaccination start!\n");
    exit(1);
  }else{
    sscanf(char_tstartvax_option,"%d",&param.tstartvax);
  }
  fprintf(stderr,"tstartvax=%d\n",param.tstartvax);

  /*time at which transmission starts*/
  if((char_tstartinf_option=get_value(argv,argc,"-tstartinf"))==NULL){
    fprintf(stderr,"Specify day of epidemic seeding!\n");
    exit(1);
  }else{
    sscanf(char_tstartinf_option,"%d",&param.tstartinf);
  }
  fprintf(stderr,"tstartinf=%d\n",param.tstartinf);

   /*target vaccination coverage (selector for file with number of doses)*/
  if((char_cov_option=get_value(argv,argc,"-cov"))==NULL){
    fprintf(stderr,"Specify vaccination coverage!\n");
    exit(1);
  }else{
    sscanf(char_cov_option,"%d",&param.cov);
  }
  fprintf(stderr,"target coverage=%d\n",param.cov);
  
  /*system capacity, in terms of first doses, i.e. half the overall capacity of the system (selector for file with number of doses)*/
  if((char_capacity_option=get_value(argv,argc,"-capacity"))==NULL){
    fprintf(stderr,"Specify capacity!\n");
    exit(1);
  }else{
    sscanf(char_capacity_option,"%d",&param.capacity);
  }
  fprintf(stderr,"system capacity=%d\n",param.capacity);

  /*omega_1: transition rate between V0 and V1*/
  if((char_capacity_option=get_value(argv,argc,"-omega_1"))==NULL){
    fprintf(stderr,"Specify omega_1!\n");
    exit(1);
  }else{
    sscanf(char_capacity_option,"%lf",&param.omega_1);
  }
  fprintf(stderr,"omega_1=%f\n",param.omega_1);
  
  /*omega_2: transition rate between V1 and V2*/
  if((char_capacity_option=get_value(argv,argc,"-omega_2"))==NULL){
    fprintf(stderr,"Specify omega_2!\n");
    exit(1);
  }else{
    sscanf(char_capacity_option,"%lf",&param.omega_2);
  }
  fprintf(stderr,"omega_2=%f\n",param.omega_2);

  /*waning_rate: transition rate between V2 and W*/
  if((char_capacity_option=get_value(argv,argc,"-waning_rate"))==NULL){
    fprintf(stderr,"Specify waning_rate!\n");
    exit(1);
  }else{
    sscanf(char_capacity_option,"%lf",&param.waning_rate);
  }
  fprintf(stderr,"waning_rate=%f\n",param.waning_rate);

  
  /*-----------------------------------------------------*/
  /* Initialize random seed and create output directory  */
  /*-----------------------------------------------------*/

  /*initialize random seed*/
  setenv("GSL_RNG_SEED",  char_seed, 1);
  gsl_rng_env_setup();
  R_GLOBAL = gsl_rng_alloc (gsl_rng_default);
  
  /*create output directories*/
  sprintf(EXP_DIR,"%s",EXP_DIR);
  sprintf(outdir,"%s/%s",EXP_DIR, char_where); 
  mkdir(outdir,S_IRWXU);
  sprintf(outdir,"%s/%s/interv_%d",EXP_DIR,  char_where,param.interv);

  mkdir(outdir,S_IRWXU);

  /*-----------------------------------------------------*/
  /*Read input files */
  /*-----------------------------------------------------*/
 
  /*read parameters file*/
  sprintf(param_file,"%s/../common_input/parameters_%s",EXP_DIR,char_where);

  if(read_param(param_file,'\t', &param)!=0){
    fprintf(stderr,"error in read_param\n");
    exit(1);
  }

  /*recovery rate from infection. 1/gamma=average generation time*/
  fprintf(stderr,"gamma=%f\n",param.gamma);
  /*number of iterations*/
  fprintf(stderr,"Nit=%d\n",param.Nit);
  /*number of days simulated*/
  fprintf(stderr,"Tmax=%d\n",param.Tmax);
  /*number of time steps within one day*/
  fprintf(stderr,"ZETA=%d\n",param.ZETA);
  /*number of vaccination compartments*/
  fprintf(stderr,"nVAX_COMPARTMENT=%d\n",NVAXCOMP);


  /*time step*/
  param.DELTAT=1/(param.ZETA*1.0);
  fprintf(stderr,"DELTAT=%f\n",param.DELTAT);
  /*transform rates in probabilities in each time step*/
  param.prob_gamma=1.-exp(-param.gamma*param.DELTAT);
  param.prob_omega_1=1.-exp(-param.omega_1*param.DELTAT);
  param.prob_omega_2=1.-exp(-param.omega_2*param.DELTAT);
  param.prob_waning=1.-exp(-param.waning_rate*param.DELTAT);
  
  /*read vector of attenuated transmission rates */
  sprintf(i_file,"%s/../common_input/betas/%s/immunity_%s/beta_%s_SUSC_%d_Reff_%s",EXP_DIR,char_where,char_iscen,char_where,param.is,char_phi);

  
  input_file=fopen(i_file,"r"); 
  if(input_file==NULL){ 
    fprintf(stderr,"ERROR: file %s not found\n",i_file); 
    exit(1);
  }

  param.beta=(double*)calloc(param.Nit,sizeof(double));  
  for(i=0;i<param.Nit;i++)
    fscanf(input_file,"%lf\n", &(param.beta[i])); 
  fclose(input_file); 
  
  /*read population by age with underlying conditions in the simulated state*/
  sprintf(age_file,"%s/../common_input/age_structure/age_structure%s_u",EXP_DIR, char_where);
  if(read_age_structure(age_file,'\t',&param)!=0){
    fprintf(stderr,"error in read_age_structure\n");
    exit(1);
  }
  
  param.N_u=(int*)calloc(param.n_ages,sizeof(int));
  for(j=0; j<param.n_ages; j++){
    param.N_u[j]=param.N[j];
  }
 
  /* read population by age without underlying conditions in the simulated state */
  sprintf(age_file,"%s/../common_input/age_structure/age_structure%s_nu",EXP_DIR, char_where);
  if(read_age_structure(age_file,'\t',&param)!=0){
    fprintf(stderr,"error in read_age_structure\n");
    exit(1);
  }
  
  param.N_nu=(int*)calloc(param.n_ages,sizeof(int));
  for(j=0; j<param.n_ages; j++){
    param.N_nu[j]=param.N[j];
  }

  /*Total population by age*/
  for(j=0; j<param.n_ages; j++){
    param.N[j]=param.N_u[j]+param.N_nu[j];
  }
  
  /* Read vector of age-specific probability of developing symptoms (same for u/nu) */
  sprintf(i_file,"%s/../common_input/psym",EXP_DIR);
  
  input_file=fopen(i_file,"r"); 
  if(input_file==NULL){ 
    fprintf(stderr,"ERROR: file %s not found\n",i_file); 
    exit(1);
  }

  param.psym=(double*)calloc(param.n_ages,sizeof(double));  
  for(i=0;i<param.n_ages;i++)
    fscanf(input_file,"%lf\n", &(param.psym[i])); 
  fclose(input_file); 

  if(param.vax_only_susc==1){
    for(i=0;i<param.n_ages;i++){
      param.psym[i]=1;
    }
    param.perc_sym_detected=1.;
    fprintf(stderr,"update pdetection=%f\n",param.perc_sym_detected);
    fprintf(stderr,"update psym...\n");    
  }

  
  /* Read matrix of age-specific susceptibility to infection */
  param.rel_sus=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.rel_sus[j]=(double *)calloc(param.Nit,sizeof(double));

  sprintf(rs_file,"%s/../common_input/relative_susceptibility_matrix_susc_%d",EXP_DIR,param.is);
  strcpy(char_rs, "");

  
  if(read_relative_susceptibility_uncertainty(rs_file,'\t',&param)!=0){
    fprintf(stderr,"error in read_relative_susceptibility\n");
    exit(1);
  }

  /* Read matrix of age-specific initial immunity to infection */
  param.prob_immune=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.prob_immune[j]=(double *)calloc(param.Nit,sizeof(double));
  
  sprintf(rs_file,"%s/../common_input/initial_immunity/%s/initial_immunity_%s",EXP_DIR,char_where,char_iscen);
  strcpy(char_rs, "");
  
  if(read_initial_immunity_uncertainty(rs_file,'\t',&param)!=0){
    fprintf(stderr,"error in read_initial_immunity\n");
    exit(1);
  }

  /*Read number of people to be vaccinated by age in each time step according to the simulated intervention */

  /*support matrix to read input_files*/
  param.vax_cov=(int **)calloc(param.Tmax*param.ZETA+1,sizeof(int *));

  /*number of doses to be administered  in each time step to individuals with underlying conditions*/
  param.vax_cov_u=(int **)calloc(param.Tmax*param.ZETA+1,sizeof(int *));
  param.vax_cov_tmp_u=(int **)calloc(param.Tmax*param.ZETA+1,sizeof(int *)); /*for checks*/

  /*number of doses to be administered  in each time step to individuals without underlying conditions*/
  param.vax_cov_nu=(int **)calloc(param.Tmax*param.ZETA+1,sizeof(int *));
  param.vax_cov_tmp_nu=(int **)calloc(param.Tmax*param.ZETA+1,sizeof(int *)); /*for checks*/

  for(i=0; i<=param.Tmax*param.ZETA; i++){
    param.vax_cov[i]=(int *)calloc(param.n_ages,sizeof(int));
    param.vax_cov_u[i]=(int *)calloc(param.n_ages,sizeof(int));
    param.vax_cov_tmp_u[i]=(int *)calloc(param.n_ages,sizeof(int));
    param.vax_cov_nu[i]=(int *)calloc(param.n_ages,sizeof(int));
    param.vax_cov_tmp_nu[i]=(int *)calloc(param.n_ages,sizeof(int));
  }
  
  /* read input doses underlying conditions */
  sprintf(cm_file,"%s/../common_input/vaccination/ndoses/%s/interv_%d/u_ndoses_%s_capacity_%d_cov_%d",EXP_DIR,char_where,param.interv,char_where,param.capacity,param.cov);
  if(read_vaccination_coverage(cm_file,'\t',&param)!=0){ 
    fprintf(stderr,"error in read_contact_matrix\n");
    exit(1);
  }

  for(i=1; i<=param.Tmax*param.ZETA; i++){
    for(j=0; j<param.n_ages; j++){
      param.vax_cov_u[i][j]=param.vax_cov_tmp_u[i][j]=param.vax_cov[i][j];
    }
  }
  
  /* read input doses without underlying conditions */
  sprintf(cm_file,"%s/../common_input/vaccination/ndoses/%s/interv_%d/nu_ndoses_%s_capacity_%d_cov_%d",EXP_DIR,char_where,param.interv,char_where,param.capacity,param.cov);
  if(read_vaccination_coverage(cm_file,'\t',&param)!=0){ 
    fprintf(stderr,"error in read_contact_matrix\n");
    exit(1);
  }

  for(i=1; i<=param.Tmax*param.ZETA; i++){
    for(j=0; j<param.n_ages; j++){
       param.vax_cov_nu[i][j]=param.vax_cov_tmp_nu[i][j]=param.vax_cov[i][j];
    }
  }
  
  /* if selector vaccine efficacy is 0, set number of doses to 0 (no vaccination)*/
  if(strcmp(char_vaxeff,"0")==0){
    for(i=0; i<=(param.Tmax*param.ZETA); i++){
      for(j=0; j<param.n_ages; j++){
	param.vax_cov_u[i][j]=0;
	param.vax_cov_tmp_u[i][j]=0;
	param.vax_cov_nu[i][j]=0;
	param.vax_cov_tmp_nu[i][j]=0;
      }
    }
    fprintf(stderr,"scenario with vaccine efficacy 0 corresponds to no vaccination\n Set number of daily doses to 0\n\n");
  }
      
  /*Allocation of variables */
  
  wasted_capacity=(int *)calloc(param.Tmax*param.ZETA+1,sizeof(int));
  wasted_capacity_u=(int *)calloc(param.Tmax*param.ZETA+1,sizeof(int));
  wasted_capacity_nu=(int *)calloc(param.Tmax*param.ZETA+1,sizeof(int));

  second_doses_over_time=(int **)calloc(param.Tmax,sizeof(int*));
  for(i=0; i<param.Tmax; i++){
     second_doses_over_time[i]=(int *)calloc(param.Nit,sizeof(int));
  }
  
  cumdet_u=(int *)calloc(param.n_ages,sizeof(int));
  cumundet_u=(int *)calloc(param.n_ages,sizeof(int));
  cumdet_nu=(int *)calloc(param.n_ages,sizeof(int));
  cumundet_nu=(int *)calloc(param.n_ages,sizeof(int));
  
  cuminf=(int *)calloc(param.n_ages,sizeof(int));
  
  /*Read age-specific vaccine efficacy in the different ramp-up stages*/
  /*Note: number of columns of vaccine efficacy should coincide with number of vaccination compartments (NVAXCOMP)*/

  param.vax_eff=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.vax_eff[j]=(double *)calloc(NVAXCOMP,sizeof(double));
  
  sprintf(rs_file,"%s/../common_input/vaccination/vaccine_efficacy/vaccine_efficacy_%s",
	  EXP_DIR,char_vaxeff);
  strcpy(char_rs, "");

  
  if(read_vaccine_efficacy(rs_file,'\t',&param)!=0){
    fprintf(stderr,"error in read_vaccine_efficacy\n");
    exit(1);
  }

  if(param.vax_prevent==1){ /*if vaccination prevents only symptoms, set vaccine efficacy to 0*/
    fprintf(stderr,"Set vaccine efficacy to 0...\n");
    for(i=0;i<param.n_ages;i++){
      for(j=0; j<NVAXCOMP; j++){
	param.vax_eff[i][j]=0.;
      }
    }
  }
  /*Allocate contact matrices*/
  param.CMhouse=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.CMhouse[j]=(double *)calloc(param.n_ages,sizeof(double));
  
  param.CMschool=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.CMschool[j]=(double *)calloc(param.n_ages,sizeof(double));
  
  param.CMworkall=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.CMworkall[j]=(double *)calloc(param.n_ages,sizeof(double));

  param.CMrandom=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.CMrandom[j]=(double *)calloc(param.n_ages,sizeof(double));
    
  param.CMtot=(double **)calloc(param.n_ages,sizeof(double *));
  for(j=0; j<param.n_ages; j++)
    param.CMtot[j]=(double *)calloc(param.n_ages,sizeof(double));
  
 
  /* Loop iterations*/
  for(i=0;i<param.Nit;i++){
    CURRENT_SIMULATION=i+1;
   
    /* For each iteration read from file the corresponding bootstapped contact matrices by setting */

    /*House*/
    sprintf(cm_file,"%s/../common_input/contact_matrices/%s/ac_h_ext_%d",EXP_DIR,char_where,CURRENT_SIMULATION);
    if(read_contact_matrix_by_setting(cm_file,'\t',&param, HOUSE_IDX)!=0){ 
      fprintf(stderr,"error in read_contact_matrix\n");
      exit(1);
    }
    /*School*/
    sprintf(cm_file,"%s/../common_input/contact_matrices/%s/ac_s_ext_%d",EXP_DIR,char_where,CURRENT_SIMULATION);
    if(read_contact_matrix_by_setting(cm_file,'\t',&param, SCHOOL_IDX)!=0){
      fprintf(stderr,"error in read_contact_matrix\n");
      exit(1);
    }
    /*Work*/
    sprintf(cm_file,"%s/../common_input/contact_matrices/%s/ac_w_ext_%d",EXP_DIR,char_where,CURRENT_SIMULATION);
    if(read_contact_matrix_by_setting(cm_file,'\t',&param, WORK_IDX)!=0){
      fprintf(stderr,"error in read_contact_matrix\n");
      exit(1);
    }
    /*Random*/
    sprintf(cm_file,"%s/../common_input/contact_matrices/%s/ac_r_ext_%d",EXP_DIR,char_where,CURRENT_SIMULATION);
    if(read_contact_matrix_by_setting(cm_file,'\t',&param, RANDOM_IDX)!=0){
      fprintf(stderr,"error in read_contact_matrix\n");
      exit(1);
    }
   
    /*Incorporate age-specific susceptibility r_a into the contact matrices*/
    for(j=0; j<param.n_ages; j++){
      for(q=0; q<param.n_ages; q++){
	param.CMhouse[j][q]*=param.rel_sus[j][i];
	param.CMschool[j][q]*=param.rel_sus[j][i];
	param.CMworkall[j][q]*=param.rel_sus[j][i];
	param.CMrandom[j][q]*=param.rel_sus[j][i];
      }
    }

        
    build_total_matrix(&param); 


    /*open output file of new daily infections by age*/
    sprintf(output_file,"%s/resA_%s_phi%s_vaxeff%s_susc%d_interv%d_iscen%s_capacity%d_cov%d_sim_%d.tsv", outdir,char_where,char_phi,char_vaxeff, param.is, param.interv,char_iscen,param.capacity,param.cov,CURRENT_SIMULATION);
    
    fpout=fopen(output_file, "w");   
    if(fpout==NULL){
      fprintf(stderr,"Error opening file cases for writing\n");
      exit(1);
    }

    /*If vaccination prevent only symptoms, the code keeps track of infections generated by the compartment of effectively vaccinated individuals (V2) - to these individuals a lower probability of developing symptoms will be applied*/
    	 
    if(param.vax_prevent==1){
      sprintf(output_file,"%s/resB_%s_phi%s_vaxeff%s_susc%d_interv%d_iscen%s_capacity%d_cov%d_sim_%d.tsv", outdir,char_where,char_phi,char_vaxeff, param.is, param.interv,char_iscen,param.capacity,param.cov,CURRENT_SIMULATION);
      fpout_vax=fopen(output_file, "w");   
      if(fpout_vax==NULL){
	fprintf(stderr,"Error opening file cases for writing\n");
	exit(1);
      }
    }
    
    
    
    /* Allocate state variables for population with underlying conditions*/
    C_u=(int*) calloc(param.n_ages, sizeof(int)); /* counter of the daily number of new infections by age */
    C_u_vax=(int*) calloc(param.n_ages, sizeof(int)); /* counter of the daily number of new infections by age among effectively vaccinated */ 

    S_u=(int*) calloc(param.n_ages, sizeof(int)); /* susceptibles by age */
    I_u=(int*) calloc(param.n_ages, sizeof(int)); /* infectious by age */
    R_u=(int*) calloc(param.n_ages, sizeof(int)); /* recovered by age */
    V0_u=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated 1st dose by age */
    V1_u=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated 2nd dose by age (not effective yet) */
    V2_u=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated 2nd dose by age (effective) */
    W_u=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated for whom immunity has waned */
  
    /* Allocate state variables for population without underlying conditions*/
    C_nu=(int*) calloc(param.n_ages, sizeof(int)); /* counter of the daily number of new infections by age */
    C_nu_vax=(int*) calloc(param.n_ages, sizeof(int)); /* counter of the daily number of new infections by age among effectively vaccinated */ 

    S_nu=(int*) calloc(param.n_ages, sizeof(int)); /* susceptibles by age */
    I_nu=(int*) calloc(param.n_ages, sizeof(int)); /* infectious by age*/
    R_nu=(int*) calloc(param.n_ages, sizeof(int)); /* recovered by age */
    V0_nu=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated 1st dose by age */
    V1_nu=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated 2nd dose by age (not effective yet)*/
    V2_nu=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated 2nd dose by age (effective) */
    W_nu=(int*) calloc(param.n_ages, sizeof(int)); /* vaccinated for whom immunity has waned */

    
    nsecond_doses_u=(int*) calloc(param.n_ages, sizeof(int)); 
    nsecond_doses_nu=(int*) calloc(param.n_ages, sizeof(int)); 

    /*Initialize variables*/
    for(j=0; j<param.n_ages; j++){
      cumundet_u[j]=0;
      cumundet_nu[j]=0;

      cumdet_u[j]=0;
      cumdet_nu[j]=0;
      cuminf[j]=0;
    }

     /*Initialize immunological status of the population*/
    for(j=0; j<param.n_ages; j++){
      /*underlying*/
      S_u[j]=round(param.N_u[j]*(1.-param.prob_immune[j][i])); /*susceptibles*/
      R_u[j]=param.N_u[j]-S_u[j]; /*removed*/
      if(strcmp(char_iscen,"0")!=0){ 
	/*compute cumulative number of undetected infections occurred before the start of vaccination*/
	/*these individuals will be eligible for vaccination*/
	cumundet_u[j]+= gsl_ran_binomial(R_GLOBAL,
				      1.-param.perc_sym_detected*param.psym[j],
				      R_u[j]); /*number of undetected infections*/
      }
      /*not underlying*/
      S_nu[j]=round(param.N_nu[j]*(1.-param.prob_immune[j][i]));/*susceptibles*/
      R_nu[j]=param.N_nu[j]-S_nu[j]; /*removed*/
      if(strcmp(char_iscen,"0")!=0){
	/*compute cumulative number of undetected infections occurred before the start of vaccination*/
	/*these individuals will be eligible for vaccination*/
	cumundet_nu[j]+= gsl_ran_binomial(R_GLOBAL,
				      1.-param.perc_sym_detected*param.psym[j],
				      R_nu[j]); /*number of undetected infections*/
      }
    }


 
    for(t=1; t<=param.Tmax*param.ZETA; t++){
      wasted_capacity[t]=0;
      wasted_capacity_u[t]=0;
      wasted_capacity_nu[t]=0;
      for(j=0; j<param.n_ages; j++){
	param.vax_cov_tmp_u[t][j]=param.vax_cov_u[t][j];
	param.vax_cov_tmp_nu[t][j]=param.vax_cov_nu[t][j];
      }
    }
    
    

    T=0;
    int first_time_step=0; /*inital infections are seeded the first time step of T=tstartinf */
    /*Start of dynamic simulations*/
    for(t=1; t<=(param.Tmax*param.ZETA); t++){
      if(T==param.tstartinf){
	/*Seed initial infections*/
	if(first_time_step==0){
	  for(j=0;j<initial_infections;j++){
	    q=gsl_rng_uniform_int(R_GLOBAL,param.n_ages);
	    while(S_nu[q]==0){
	      q=gsl_rng_uniform_int(R_GLOBAL,param.n_ages);
	    }
	    S_nu[q]--;
	    I_nu[q]++;
	    C_nu[q]++;
	  
	    if(S_nu[q]<0){
	      fprintf(stderr,"Error!\n");
	      exit(1);
	    }
	  }
	}
	first_time_step=1;
      }

      /*Allocate counters new events*/
      nI_u=(int *) calloc(param.n_ages, sizeof(int));
      nIdet_u=(int *) calloc(param.n_ages, sizeof(int));
      nI_0_u=(int *) calloc(param.n_ages, sizeof(int));
      nI_1_u=(int *) calloc(param.n_ages, sizeof(int));
      nI_2_u=(int *) calloc(param.n_ages, sizeof(int));
      nI_3_u=(int *) calloc(param.n_ages, sizeof(int));
    
      nR_u=(int *) calloc(param.n_ages, sizeof(int));
      nV0_u=(int *) calloc(param.n_ages, sizeof(int));
      nV1_u=(int *) calloc(param.n_ages, sizeof(int));
      nV2_u=(int *) calloc(param.n_ages, sizeof(int));
      nW_u=(int *) calloc(param.n_ages, sizeof(int));

      nI_nu=(int *) calloc(param.n_ages, sizeof(int));
      nIdet_nu=(int *) calloc(param.n_ages, sizeof(int));
      nI_0_nu=(int *) calloc(param.n_ages, sizeof(int));
      nI_1_nu=(int *) calloc(param.n_ages, sizeof(int));
      nI_2_nu=(int *) calloc(param.n_ages, sizeof(int));
      nI_3_nu=(int *) calloc(param.n_ages, sizeof(int));
    
      nR_nu=(int *) calloc(param.n_ages, sizeof(int));
      nV0_nu=(int *) calloc(param.n_ages, sizeof(int));
      nV1_nu=(int *) calloc(param.n_ages, sizeof(int));
      nV2_nu=(int *) calloc(param.n_ages, sizeof(int));
      nW_nu=(int *) calloc(param.n_ages, sizeof(int));
     
      
      /*Allocate force of infection by age (same for those with/without u.c.)*/
      foi=(double*) calloc(param.n_ages, sizeof(double));
      /*Initialize age-dependent force of infection*/ 
      for(j=0; j<param.n_ages; j++){
	foi[j]=0.0;
      }
	
      /*Compute age-dependent force of infection*/ 
      for(j=0; j<param.n_ages; j++){
	for(q=0; q<param.n_ages; q++){
	  foi[j]+=
	    param.CMtot[j][q]*param.beta[i]*(I_u[q]+I_nu[q])/ (double) (param.N_u[q]+param.N_nu[q]);
	}
      }

  
      /*Epidemiological transitions*/
      for(j=0; j<param.n_ages; j++){
	
	/*1) vaccination*/
	if(T>=param.tstartvax){ 

	  /* 1A) Vaccination of people with underlying conditions*/
	  if(param.vax_cov_u[t-(param.tstartvax*param.ZETA)][j]>0){
	    
	    if(S_u[j]>0 || cumundet_u[j]>0){
	      /*distribute first doses among susceptibles and infected/removed undetected (and not vaccinated yet)*/
	      /*first doses to be "theoretically" administered to susceptibles*/
	      nvax_u=param.vax_cov_u[t-(param.tstartvax*param.ZETA)][j]*(S_u[j]/(1.*(S_u[j]+cumundet_u[j])));
	      /*first doses to be "theoretically" administered to undetected (wasted, because they do not prevent from infection, these individuals have already developed SARS-CoV-2 infeciton) */
	      nvax_wasted_u=param.vax_cov_u[t-(param.tstartvax*param.ZETA)][j]-nvax_u;
	    }else{
	      nvax_u=0;
	      nvax_wasted_u=0;	      
	    }

	    /*if number of doses is greater than available susceptibles, administer all possible doses (rest of daily capacity is wasted)*/
	    if(nvax_u>S_u[j]){
	      wasted_capacity_u[t-(param.tstartvax*param.ZETA)]+=(nvax_u-S_u[j]);
	      nvax_u=S_u[j];
	    }

	    /*if number of doses is greater than available undetected, administer all possible doses (rest of daily capacity is wasted)*/
	    if(nvax_wasted_u>cumundet_u[j]){
	      wasted_capacity_u[t-(param.tstartvax*param.ZETA)]+=(nvax_wasted_u-cumundet_u[j]);
	      nvax_wasted_u=cumundet_u[j];
	    }
	    
	    cumundet_u[j]-=nvax_wasted_u; /*subtract from counter undtected that received vaccination to "avoid" re-vaccination*/

	    param.vax_cov_tmp_u[t-(param.tstartvax*param.ZETA)][j]-=(nvax_u+nvax_wasted_u);
	    if(param.vax_cov_tmp_u[t-(param.tstartvax*param.ZETA)][j]<0){
	      fprintf(stderr,"t=%d\tERROR!!",t);
	      exit(1);
	    }
	  }else{
	    nvax_u=0;
	    nvax_wasted_u=0;	    
	  }

	  /* 1B) Vaccination of people with underlying conditions*/
	  if(param.vax_cov_nu[t-(param.tstartvax*param.ZETA)][j]>0){

	    if(S_nu[j]>0 || cumundet_nu[j]>0){
	      nvax_nu=param.vax_cov_nu[t-(param.tstartvax*param.ZETA)][j]*(S_nu[j]/(1.*(S_nu[j]+cumundet_nu[j])));
	      
	      nvax_wasted_nu=param.vax_cov_nu[t-(param.tstartvax*param.ZETA)][j]-nvax_nu;
	    }else{
	      nvax_nu=0;
	      nvax_wasted_nu=0;	      
	    }

	    /*if number of doses is greater than available susceptibles, administer all possible doses (rest of daily capacity is wasted)*/
	    if(nvax_nu>S_nu[j]){
	      wasted_capacity_nu[t-(param.tstartvax*param.ZETA)]+=(nvax_nu-S_nu[j]);
	      nvax_nu=S_nu[j];
	    }

	    /*if number of doses is greater than available undetected, administer all possible doses (rest of daily capacity is wasted)*/
	    if(nvax_wasted_nu>cumundet_nu[j]){
	      wasted_capacity_nu[t-(param.tstartvax*param.ZETA)]+=(nvax_wasted_nu-cumundet_nu[j]);
	      nvax_wasted_nu=cumundet_nu[j];
	    }
	    
	    cumundet_nu[j]-=nvax_wasted_nu; /*subtract from counter undtected that received vaccination to "avoid" re-vaccination*/
	    param.vax_cov_tmp_nu[t-(param.tstartvax*param.ZETA)][j]-=(nvax_nu+nvax_wasted_nu);
	    
	    if(param.vax_cov_tmp_nu[t-(param.tstartvax*param.ZETA)][j]<0){
	      fprintf(stderr,"t=%d\tERROR!!",t);
	      exit(1);
	    }
	  }else{
	    nvax_nu=0;
	    nvax_wasted_nu=0;
	  }

	  /*New transitions between S and V0*/
	  nV0_u[j]=nvax_u;
	  nV0_nu[j]=nvax_nu;
	 
	  
	  /*progression between the different ramp-up stages*/
	  nV1_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_omega_1, V0_u[j]);
	  nV1_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_omega_1, V0_nu[j]);
	  
	  nV2_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_omega_2 , V1_u[j]);
	  nV2_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_omega_2 , V1_nu[j]);

	  nW_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_waning , V2_u[j]);
	  nW_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_waning , V2_nu[j]);

	  nsecond_doses_u[j]+=(nV1_u[j]);
	  nsecond_doses_nu[j]+=(nV1_nu[j]);

	}

	if(T>=param.tstartinf){ /*start transmission*/
	  /* force of infection to which are subject S and W */
	  param.prob_foi=1.-exp(-foi[j]*param.DELTAT);
	  /* force of infection to which are subject V0 */
	  param.prob_foi_0=1.-exp(-(1-param.vax_eff[j][0])*foi[j]*param.DELTAT);
	  /* force of infection to which are subject V1 */
	  param.prob_foi_1=1.-exp(-(1-param.vax_eff[j][1])*foi[j]*param.DELTAT);
	  /* force of infection to which are subject V2 */
	  param.prob_foi_2=1.-exp(-(1-param.vax_eff[j][2])*foi[j]*param.DELTAT); 
	}else{
	  /*no new infections*/
	  param.prob_foi=0.; 
	  param.prob_foi_0=0.; 
	  param.prob_foi_1=0.; 
	  param.prob_foi_2=0.;
	}
	
	/*2) new infections*/
	/* - among unvaccinated with or without u.c. */
	nI_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi, S_u[j]-nV0_u[j]); /*among unvaccinated with u.c. */
	nI_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi, S_nu[j]-nV0_nu[j]); /*among unvaccinated with u.c. */
	
	
	nIdet_u[j]=gsl_ran_binomial(R_GLOBAL, param.perc_sym_detected*param.psym[j], nI_u[j]); /*number of detected infections*/
	nIdet_nu[j]=gsl_ran_binomial(R_GLOBAL, param.perc_sym_detected*param.psym[j], nI_nu[j]); /*number of detected infections*/

	cuminf[j]+=(nI_u[j]+nI_nu[j]); /*used only for checks*/
	
	cumdet_u[j]+=nIdet_u[j]; /*update counter detected infections (these  will be excluded from vaccination)*/
	cumdet_nu[j]+=nIdet_nu[j]; /*update counter detected infections (these  will be excluded from vaccination)*/

	/*update counters undetected infections*/
	cumundet_u[j]+=(nI_u[j]-nIdet_u[j]);
	cumundet_nu[j]+=(nI_nu[j]-nIdet_nu[j]);

	/* - among vaccinated with or without u.c. */
	
	nI_0_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi_0, V0_u[j]+nV0_u[j]-nV1_u[j]); 
	nI_1_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi_1, V1_u[j]+nV1_u[j]-nV2_u[j]);
	nI_2_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi_2, V2_u[j]+nV2_u[j]-nW_u[j]);
	nI_3_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi, W_u[j]+nW_u[j]);
	
	nI_0_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi_0, V0_nu[j]+nV0_nu[j]-nV1_nu[j]); 
	nI_1_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi_1, V1_nu[j]+nV1_nu[j]-nV2_nu[j]);
	nI_2_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi_2, V2_nu[j]+nV2_nu[j]-nW_nu[j]);
	nI_3_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_foi, W_nu[j]+nW_nu[j]);


	/*3) recoveries */
	nR_u[j]=gsl_ran_binomial(R_GLOBAL, param.prob_gamma, I_u[j]);
	nR_nu[j]=gsl_ran_binomial(R_GLOBAL, param.prob_gamma, I_nu[j]);


	/*Update state variables*/
	S_u[j]=S_u[j]-nI_u[j]-nV0_u[j];
	V0_u[j]=V0_u[j]+nV0_u[j]-nI_0_u[j]-nV1_u[j];	
	V1_u[j]=V1_u[j]+nV1_u[j]-nI_1_u[j]-nV2_u[j];
	V2_u[j]=V2_u[j]+nV2_u[j]-nI_2_u[j]-nW_u[j];
	W_u[j]=W_u[j]+nW_u[j]-nI_3_u[j];
	I_u[j]=I_u[j]+nI_u[j]+nI_0_u[j]+nI_1_u[j]+nI_2_u[j]+nI_3_u[j]-nR_u[j];
	R_u[j]=R_u[j]+nR_u[j];

	S_nu[j]=S_nu[j]-nI_nu[j]-nV0_nu[j];
	V0_nu[j]=V0_nu[j]+nV0_nu[j]-nI_0_nu[j]-nV1_nu[j];	
	V1_nu[j]=V1_nu[j]+nV1_nu[j]-nI_1_nu[j]-nV2_nu[j];
	V2_nu[j]=V2_nu[j]+nV2_nu[j]-nI_2_nu[j]-nW_nu[j];
	W_nu[j]=W_nu[j]+nW_nu[j]-nI_3_nu[j];
	I_nu[j]=I_nu[j]+nI_nu[j]+nI_0_nu[j]+nI_1_nu[j]+nI_2_nu[j]+nI_3_nu[j]-nR_nu[j];
	R_nu[j]=R_nu[j]+nR_nu[j];
	
	/*update counter new infections by age*/ 
	C_u[j]+=(nI_u[j]+nI_0_u[j]+nI_1_u[j]+nI_2_u[j]+nI_3_u[j]); 
	C_nu[j]+=(nI_nu[j]+nI_0_nu[j]+nI_1_nu[j]+nI_2_nu[j]+nI_3_nu[j]); 

	/*update counter new infections by age among effectively vaccinated*/ 
	C_u_vax[j]+=(nI_2_u[j]); 
	C_nu_vax[j]+=(nI_2_nu[j]); 

	
	/*Check 1*/
	if(S_u[j]<0 || I_u[j]<0 || R_u[j]<0 || V0_u[j]<0 || V1_u[j]<0 || V2_u[j]<0 || W_u[j]<0 ||
	   S_nu[j]<0 || I_nu[j]<0 || R_nu[j]<0 || V0_nu[j]<0 || V1_nu[j]<0 || V2_nu[j]<0 || W_nu[j]<0){
	  fprintf(stderr,"Negative population in S/I/J/K/V0/V1/V2!!\n");
	  fprintf(stderr,"t=%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		  t,S_u[j],I_u[j],R_u[j],V0_u[j],V1_u[j],V2_u[j],W_u[j],
		  S_nu[j],I_nu[j],R_nu[j],V0_nu[j],V1_nu[j],V2_nu[j],W_nu[j]);
	  exit(1);
	}

	/*Check 2*/
	if(strcmp(char_vaxeff,"0")==0 && (V0_u[j]>0 || V0_nu[j]>0 )){
	  fprintf(stderr,"ERROR!!!");
	  exit(1);
	}

	/*Check 3*/
	if((S_u[j]+I_u[j]+R_u[j]+V0_u[j]+V1_u[j]+V2_u[j]+W_u[j])!=param.N_u[j]){
	  fprintf(stderr,"Population changes over time!\n");
	  fprintf(stderr,"t=%d\tage=%d\tS=%d\tI=%d\tR=%d\tV0=%d\tV1=%d\tV2=%d\tW=%d\tS+I+R+V0+V1+V2+W=%d\tN=%d\n",
		  t,j,
		  S_u[j],I_u[j],R_u[j],V0_u[j],V1_u[j],V2_u[j],W_u[j],
		  S_u[j]+I_u[j]+R_u[j]+V0_u[j]+V1_u[j]+V2_u[j]+W_u[j],
		  param.N_u[j]);
	  exit(1);
	}
	
	/*Check 4*/
	if((S_nu[j]+I_nu[j]+R_nu[j]+V0_nu[j]+V1_nu[j]+V2_nu[j]+W_nu[j])!=param.N_nu[j]){
	  fprintf(stderr,"Population changes over time!\n");
	  fprintf(stderr,"t=%d\tage=%d\tS=%d\tI=%d\tR=%d\tV0=%d\tV1=%d\tV2=%d\tW=%d\tS+I+R+V0+V1+V2+W=%d\tN=%d\n",
		  t,j,
		  S_nu[j],I_nu[j],R_nu[j],V0_nu[j],V1_nu[j],V2_nu[j],W_nu[j],
		  S_nu[j]+I_nu[j]+R_nu[j]+V0_nu[j]+V1_nu[j]+V2_nu[j]+W_nu[j],
		  param.N_nu[j]);
	  exit(1);
	}
	
	
      }

      Itot=0;
      for(j=0; j<param.n_ages; j++){
	Itot+=(I_u[j]+I_nu[j]); 	
      }


      /*At each day, print on output files*/
      if(t%param.ZETA==0){
	for(j=0; j<param.n_ages; j++){
	  fprintf(fpout,"%d\t",C_u[j]);
	  if(param.vax_prevent==1)
	    fprintf(fpout_vax,"%d\t",C_u_vax[j]);
	  
	  second_doses_over_time[T][i]+=(nsecond_doses_u[j]+nsecond_doses_nu[j]);
	 
	}
	for(j=0; j<param.n_ages; j++){
	  fprintf(fpout,"%d\t",C_nu[j]);
	  if(param.vax_prevent==1)
	    fprintf(fpout_vax,"%d\t",C_nu_vax[j]);
	}

	fprintf(fpout, "\n");

	if(param.vax_prevent==1)
	  fprintf(fpout_vax, "\n");
	
      
	for(j=0; j<param.n_ages; j++){	  
	  C_u[j]=C_nu[j]=C_u_vax[j]=C_nu_vax[j]=0;  // reset daily counter
	  nsecond_doses_u[j]=0;  // reset daily counter
	  nsecond_doses_nu[j]=0;  // reset daily counter
	}
      }

      /*update counter days (every param.ZETA time steps)*/
      if(t%param.ZETA==0){
	T++;
       }
 
      free(foi);
      free(nI_u);
      free(nR_u);
      free(nI_0_u);
      free(nI_1_u);
      free(nI_2_u);
      free(nI_3_u);
      free(nV0_u);
      free(nV1_u);
      free(nV2_u);
      free(nW_u);

      free(nI_nu);
      free(nR_nu);
      free(nI_0_nu);
      free(nI_1_nu);
      free(nI_2_nu);
      free(nI_3_nu);
      free(nV0_nu);
      free(nV1_nu);
      free(nV2_nu);
      free(nW_nu);

      
    } /*end loop over time*/

    for(t=1; t<=(param.Tmax*param.ZETA); t++){
      wasted_capacity[t]=wasted_capacity_u[t]+wasted_capacity_nu[t];
    }


    free(nsecond_doses_u);
    free(nsecond_doses_nu);

    
    free(C_u);
    free(C_u_vax);

    free(S_u);
    free(I_u);
    free(R_u);
    free(V0_u);
    free(V1_u);
    free(V2_u);
    free(W_u);

    free(C_nu);
    free(C_nu_vax);

    free(S_nu);
    free(I_nu);
    free(R_nu);
    free(V0_nu);
    free(V1_nu);
    free(V2_nu);
    free(W_nu);

    
    fflush(fpout);
    fflush(fpout);
    fclose(fpout);

    if(param.vax_prevent==1){
      fflush(fpout_vax);
      fflush(fpout_vax);
      fclose(fpout_vax);
    }
 
  } /*end loop over number of iterations (param.Nit)*/

 
  
  
  /*free variables allocated in read_age_structure.c*/
  free(param.age_lim_lo);
  free(param.age_lim_hi);
  free(param.N);
  free(param.N_u);
  free(param.N_nu);


  /*free variables allocated in contact_matrices.c*/
  
  for(j=0; j<param.n_ages; j++){
    free(param.CMhouse[j]);
    free(param.CMschool[j]);
    free(param.CMworkall[j]);
    free(param.CMrandom[j]);
    free(param.CMtot[j]);
  }
  free(param.CMhouse);
  free(param.CMschool);
  free(param.CMworkall);
  free(param.CMrandom);
  free(param.CMtot);
  
  free(param.beta);
  free(param.psym);
  free(wasted_capacity);
  free(wasted_capacity_u);
  free(wasted_capacity_nu);
 
  
  
  for(j=0; j<param.n_ages; j++){
    free(param.rel_sus[j]);
    free(param.prob_immune[j]);
  } 
  free(param.rel_sus); 
  free(param.prob_immune); 

  
  for(j=0; j<param.n_ages; j++){
    free(param.vax_eff[j]);
  } 
  free(param.vax_eff); 


  for(j=0; j<(param.Tmax*param.ZETA); j++){
    free(param.vax_cov[j]);
    free(param.vax_cov_u[j]);
    free(param.vax_cov_nu[j]);
    free(param.vax_cov_tmp_u[j]);
    free(param.vax_cov_tmp_nu[j]);
  }
  free(param.vax_cov);
  free(param.vax_cov_u);
  free(param.vax_cov_nu);
  free(param.vax_cov_tmp_u);
  free(param.vax_cov_tmp_nu);

  for(t=0; t<param.Tmax; t++){
    free(second_doses_over_time[t]);
  }
  free(second_doses_over_time);
  
  
 
  free(cumdet_u);
  free(cumundet_u);
  free(cumdet_nu);
  free(cumundet_nu);

  free(cuminf);
  
  return 0;  
}
