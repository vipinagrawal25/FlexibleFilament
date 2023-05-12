#ifndef FILE_MESH_SEEN
#define FILE_MESH_SEEN
/*--------------------------------------------------*/
struct MESH{
  int Np=0;
  int pp=0;
  int ndim=0;
  double *y;	// address to store position
  // double *vel;	// useful for the dynamics code
  
};
/*--------------------------------------------------*/
void set_param(MESH *mesh);
#endif /* !MAPDYN_SEEN */