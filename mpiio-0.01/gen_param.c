
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char* argv[]){

  int i;
  int number_of_parameters;
  int *nsize;
  FILE *fptr;
  fptr = fopen("asci.param", "r");
  fscanf(fptr, "%d", &number_of_parameters);
  nsize = (int *)malloc(number_of_parameters * sizeof(int));

  for (i=0;i<number_of_parameters; ++i)
    fscanf(fptr, "%d", &nsize[i]);

  fclose(fptr);

  fptr = fopen("bin.param", "w");
  fwrite(&number_of_parameters,sizeof(int),1,fptr);
    for(i=0;i<number_of_parameters;++i)
          fwrite(&nsize[i],sizeof(int),1,fptr);

  fclose(fptr);

  return 0;
}
