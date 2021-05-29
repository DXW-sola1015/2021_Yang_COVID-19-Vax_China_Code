#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "header.h"
#include <math.h>

extern int HOUSE_IDX;
extern int SCHOOL_IDX;
extern int WORK_IDX;
extern int RANDOM_IDX;

 /*Read contact matrices computed from contact diaries (POLYMOD)*/
int read_contact_matrix_by_setting(char file[],char sep,Param *param,
				   int setting)
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

    for(i=0; i<param->n_ages; i++){
      if(setting==HOUSE_IDX){
	sscanf(line,"%lf", &(param->CMhouse[n_rows][i]));
      } else if(setting==SCHOOL_IDX){ 
	sscanf(line,"%lf", &(param->CMschool[n_rows][i]));
      }else if(setting==WORK_IDX){ 
	sscanf(line,"%lf", &(param->CMworkall[n_rows][i]));
      }else if(setting==RANDOM_IDX){ 
	sscanf(line,"%lf", &(param->CMrandom[n_rows][i]));
      }
      
      line = (char *)strchr(line, sep);
      line++;
    }
    n_rows++;

  }

  return 2;
}


int build_total_matrix(Param *param)
{
  int j,q;
 
  
 /*Compute matrix of overall contacts*/
 for(j=0; j<param->n_ages; j++){
    for(q=0; q<param->n_ages; q++){
      param->CMtot[j][q]=param->CMhouse[j][q]+
	param->CMschool[j][q]+
	param->CMworkall[j][q]+
	param->CMrandom[j][q];
    }
 }
 
 
  return 2;
}
