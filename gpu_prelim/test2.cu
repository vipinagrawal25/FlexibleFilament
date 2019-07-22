#include <stdio.h>
#include <string.h>
/* testing function pointers */
void rnkt4( void ) {
  printf( "rnkt4\n") ;
}
void euler( void ) {
  printf( "euler\n") ;
}

void evolve( void (*algo)( void )  ){
  algo();
}
void my_print(char *A){
  printf( " %s \n", A); 
}
/*------------------------*/
int main( void ){
  char *A = "euler";
  printf ( " %s \n", A);
  void (*ALGO)( );
  printf (" %d \n", strcmp( A, "euler" ) );
  if ( strcmp( A , "euler") == 0 )
    {
      ALGO = &euler;
    } else if ( strcmp( A, "rnkt4" ) == 0 )
    {
      ALGO = &rnkt4;
    } else
    {
      printf( " algorithm\t%s\t not coded \n", A);
      exit(1);
    }
  ALGO( );
  my_print( A);
}
