#include <iostream>
//
#include "model.h"
/********************************************/
using namespace std;
/* ----------------------------------------*/
int main(){
  MESH mesh;
  set_param(&mesh);
  cout << mesh.Np << endl;
  cout << mesh.ndim << endl;
}
/********************************************/