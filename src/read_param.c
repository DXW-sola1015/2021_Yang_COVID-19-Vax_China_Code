#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "header.h"
#include <math.h>

 /*Read parameters file*/
int read_param(char file[],char sep, Param *param){

  char *line;
  int out_get_line=2;
  FILE *fp;
  char string_who[100];
  char string_value[100];

  if(!(fp = fopen(file, "r"))){
    fprintf(stderr,"read_param: error opening file %s for reading\n",file);
    return 1;
  }

  while(out_get_line>=2){
    out_get_line=get_line(&line,fp);
    if(out_get_line<3){
      switch(out_get_line){
      case 2:
	fprintf(stderr,"read_param: line of file %s does not end in newline\n",file);
	break;
      case 1:
	fprintf(stderr,"read_param: file %s contains an empty line\n",file);
	return 1;
	break;
      case 0:
      	fclose(fp);
      	return 0;
	break;
      case -1:
	fprintf(stderr,"read_param: get_line error on file %s\n",
		file);
	return 1;
      default:
	fprintf(stderr,"read_param: unrecognized exit status of get_line on file %s\n",file);
	return 1;
	break;
      }
    }

    sscanf(line,"%s", string_who);
    line = (char *)strchr(line, sep);
    line++;
    sscanf(line,"%s", string_value);
    /*recovery rate for each stage of infection*/
    if(strcmp(string_who,"gamma")==0)
      param->gamma=atof(string_value);
    /*number of iterations*/
    if(strcmp(string_who,"Nit")==0)
      param->Nit=atoi(string_value);
    /* number of simulated days */
    if(strcmp(string_who,"Tmax")==0)
      param->Tmax=atoi(string_value);
    /* time step */
    if(strcmp(string_who,"ZETA")==0)
      param->ZETA=atoi(string_value);
  }
  return 2;
}

