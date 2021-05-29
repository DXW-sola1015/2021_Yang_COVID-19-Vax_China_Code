#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "header.h"

/*Read population by age in a given region*/
int read_age_structure(char file[],char sep,Param *param)
{
  char *line;
  int out_get_line=3;
  FILE *fp;
  int n_rows;
  int tmp;

  if(!(fp = fopen(file, "r"))){
    fprintf(stderr,"ERROR: file %s not found\n",file);
    return 1;
  }

  n_rows=0;

  param->pop=0;
  param->age_lim_lo=(int*)calloc(1,sizeof(int));
  param->age_lim_hi=(int*)calloc(1,sizeof(int));
  param->N=(int*)calloc(1,sizeof(int));
  
  while(out_get_line>=2){
    out_get_line=get_line(&line,fp);
    if(out_get_line<3){
      switch(out_get_line){
      case 2:
	fprintf(stderr,"read_age_structure: line %d of file %s does not end in newline\n",n_rows,file);
	break;
      case 1:
	fprintf(stderr,"read_age_structure: file %s contains an empty line\n",file);
	return 1;
	break;
      case 0:
	param->n_ages = n_rows;
	fclose(fp);
	return 0;
	break;
      case -1:
	fprintf(stderr,"read_age_structure: get_line error on file %s\n",
		file);
	return 1;
      default:
	fprintf(stderr,"read_age_structure: unrecognized exit status of get_line on file %s\n",file);
	return 1;
	break;
      }
    }

    sscanf(line,"%d", &(param->age_lim_lo[n_rows]));
    line = (char *)strchr(line, sep);
    line++;
    
    sscanf(line,"%d", &(param->age_lim_hi[n_rows]));
    line = (char *)strchr(line, sep);
    line++;
    
    sscanf(line,"%d", &tmp);
    if(tmp<=0){ 
      fprintf(stderr,"ERROR: population in all age groups must be greater than 0\n");
      exit(1);
    }
    param->N[n_rows]=tmp;
    param->pop+=param->N[n_rows];
    line = (char *)strchr(line, sep);
    line++;

    n_rows++;
    param->age_lim_lo=(int*)realloc(param->age_lim_lo,(n_rows+1)*sizeof(int));
    param->age_lim_hi=(int*)realloc(param->age_lim_hi,(n_rows+1)*sizeof(int));
    param->N=(int*)realloc(param->N,(n_rows+1)*sizeof(int));

  }

  return 2;
}
