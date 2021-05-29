#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "header.h"
#include <math.h>

extern int NVAXCOMP;

/* read age-specific immunity to infection */
int read_initial_immunity_uncertainty(char file[],char sep,Param *param)
{
  char *line;
  int out_get_line=3;
  FILE *fp;
  int n_rows;
  int i;

   if(!(fp = fopen(file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",file);
    return 1;
  }

  n_rows=0;
  
  while(out_get_line>=2){
    out_get_line=get_line(&line,fp);
    if(out_get_line<3){
      switch(out_get_line){
      case 2:
	fprintf(stderr,"read_contact_matrix: line %d of file %s does not end in newline\n",n_rows,file);
	break;
      case 1:
	fprintf(stderr,"read_contact_matrix: file %s contains an empty line\n",file);
	return 1;
	break;
      case 0:
	fclose(fp);
	return 0;
	break;
      case -1:
	fprintf(stderr,"read_contact_matrix: get_line error on file %s\n",
		file);
	return 1;
      default:
	fprintf(stderr,"read_contact_matrix: unrecognized exit status of get_line on file %s\n",file);
	return 1;
	break;
      }
    }

    for(i=0; i<param->Nit; i++){
      sscanf(line,"%lf", &(param->prob_immune[n_rows][i]));
      line = (char *)strchr(line, sep);
      line++;
    }

    n_rows++;
    
  }
  
  return 2;
}


/* read age-specific susceptibility to infection for all iterations*/
int read_relative_susceptibility_uncertainty(char file[],char sep,Param *param)
{
  char *line;
  int out_get_line=3;
  FILE *fp;
  int n_rows;
  int i;

   if(!(fp = fopen(file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",file);
    return 1;
  }

  n_rows=0;
  

  while(out_get_line>=2){
    out_get_line=get_line(&line,fp);
    if(out_get_line<3){
      switch(out_get_line){
      case 2:
	fprintf(stderr,"read_contact_matrix: line %d of file %s does not end in newline\n",n_rows,file);
	break;
      case 1:
	fprintf(stderr,"read_contact_matrix: file %s contains an empty line\n",file);
	return 1;
	break;
      case 0:
	fclose(fp);
	return 0;
	break;
      case -1:
	fprintf(stderr,"read_contact_matrix: get_line error on file %s\n",
		file);
	return 1;
      default:
	fprintf(stderr,"read_contact_matrix: unrecognized exit status of get_line on file %s\n",file);
	return 1;
	break;
      }
    }

    for(i=0; i<param->Nit; i++){
      sscanf(line,"%lf", &(param->rel_sus[n_rows][i]));
      line = (char *)strchr(line, sep);
      line++;
    }
    n_rows++;
    
  }
  
  return 2;
}

/* read vaccine efficacy by age*/
int read_vaccine_efficacy(char file[],char sep,Param *param)
{
  char *line;
  int out_get_line=3;
  FILE *fp;
  int n_rows;
  int i;

   if(!(fp = fopen(file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",file);
    return 1;
  }

  n_rows=0;
  

  while(out_get_line>=2){
    out_get_line=get_line(&line,fp);
    if(out_get_line<3){
      switch(out_get_line){
      case 2:
	fprintf(stderr,"read_vaccine_efficacy: line %d of file %s does not end in newline\n",n_rows,file);
	break;
      case 1:
	fprintf(stderr,"read_vaccine_efficacy: file %s contains an empty line\n",file);
	return 1;
	break;
      case 0:
	fclose(fp);
	return 0;
	break;
      case -1:
	fprintf(stderr,"read_vaccine_efficacy: get_line error on file %s\n",
		file);
	return 1;
      default:
	fprintf(stderr,"read_vaccine_efficacy: unrecognized exit status of get_line on file %s\n",file);
	return 1;
	break;
      }
    }

    for(i=0; i<NVAXCOMP; i++){
      sscanf(line,"%lf", &(param->vax_eff[n_rows][i]));
      line = (char *)strchr(line, sep);
      line++;
    }
    n_rows++;
    
  }
  
  return 2;
}



/*Read vaccination coverage by age and over time*/
int read_vaccination_coverage(char file[],char sep,Param *param)
{
  char *line;
  int out_get_line=3;
  FILE *fp;
  int n_rows;
  int i;

  if(!(fp = fopen(file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",file);
    return 1;
  }

  n_rows=0;
  
  
  while(out_get_line>=2){
    out_get_line=get_line(&line,fp);
    if(out_get_line<3){
      switch(out_get_line){
      case 2:
	fprintf(stderr,"read_vaccination_coverage: line %d of file %s does not end in newline\n",n_rows,file);
	break;
      case 1:
	fprintf(stderr,"read_vaccination_coverage: file %s contains an empty line\n",file);
	return 1;
	break;
      case 0:
	fclose(fp);
	return 0;
	break;
      case -1:
	fprintf(stderr,"read_vaccination_coverage: get_line error on file %s\n",
		file);
	return 1;
      default:
	fprintf(stderr,"read_vaccination_coverage: unrecognized exit status of get_line on file %s\n",file);
	return 1;
	break;
      }
    }

    for(i=0; i<param->n_ages; i++){
	sscanf(line,"%d", &(param->vax_cov[n_rows+1][i]));
      
      line = (char *)strchr(line, sep);
      line++;
    }
    n_rows++;
    
  }
  
  return 2;
}
