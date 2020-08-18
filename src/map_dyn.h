#ifndef FILE_MAP_DYN_SEEN
#define FILE_MAP_DYN_SEEN
// #include "main.cu"
struct MV{
  double time;    // Current time of the simulation
  double dt;      // Current dt for the adaptive scheme
  int period;      // Number of iteration of the map.
  bool istab;		// Do you want to calculate the stability of the orbit?
  bool iorbit;
  bool irel_orb;
  // int ndiag;
  // int ldiag;
};
const int size_MV = 2*sizeof( double ) + sizeof( int ) + 3*sizeof(bool);
/*--------------------------------------------------*/
#endif /* !EVOLVE_SEEN */