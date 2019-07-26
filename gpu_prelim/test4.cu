#include<stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <iostream>
using namespace std;
#include "model.h"
//#include "chain.h"
#include "cuda.h"
#define   nn 10
/* ----------------------------------------------------- */
void scpy( char to[], char from[]){
  int i=0;
  while ( (to[i] = from[i]) != '\0')
    i = i+1;
}
int main( void ){
  struct CRASH *bug;
  bug = (CRASH *)malloc( size_CRASH );
  (*bug).lstop = 1;
  scpy((*bug).message, "t" ) ;
  printf( "%s\n", (*bug).message );
}
/*---------------------------------*/


