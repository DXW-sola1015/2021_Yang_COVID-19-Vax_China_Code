/*verbosity level*/
int VERBOSE;

/*inizialization of GSL*/
gsl_rng *R_GLOBAL;

/*exp dir*/
char *EXP_DIR=NULL;

/*current simulation*/
int CURRENT_SIMULATION;

/*number of vaccination compartment*/
int NVAXCOMP=3;

/*index contact matrices*/
int HOUSE_IDX=1;
int SCHOOL_IDX=2;
int WORK_IDX=3;
int RANDOM_IDX=4;
