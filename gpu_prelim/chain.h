#ifndef FILE_CHAIN_SEEN
#define FILE_CHAIN_SEEN
/*--------------------------------------------------*/
/* We solve for a chain, a one dimensional system. 
Each node may be connected to every other node. 
The number of nodes is NN. 
The number of degrees of freedom at each node is pp */ 
#define NN 16
#define pp  3
#define ndim NN*pp
void alloc_chain( double **PSI, double **psi);
void free_chain( double **PSI, double **psi );
void H2D(double psi[], double PSI[], int Nsize );
void D2H(double PSI[], double psi[], int Nsize );
void iniconf(  double PSI[], double psi[]);
#endif /* !CHAIN_SEEN */
