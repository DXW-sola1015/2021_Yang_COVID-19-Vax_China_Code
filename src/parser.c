#include <stdio.h>
#include <string.h>
#include <stdlib.h>


char *get_value(char *cl[],int ncl,char opt[]){ //GG -- questa funzione puo' essere accorpata nel file utilities.c
  int i;
  char *out=NULL;
  
  for(i=2;i<ncl-1;i++)
    if(strcmp(cl[i],opt)==0)
      out=cl[i+1];
  
  return out;
}
