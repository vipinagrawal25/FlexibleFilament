/* This file is just a copy of solve.cpp file but parallelized using MPI.*/
#include <iostream>
#include <time.h>
#include <math.h>
#include <string>
#include <vector>
#include "ode.h"
#include "model.h"
#include "input.h"
#include <sys/stat.h>
#include "constant.h"
#include "io.h"
#include "misc.h"
#include <fstream>
#include <unistd.h>
#include <memory.h>
#include <mpi.h>
/********************************************/
using namespace std;
/* ----------------------------------------*/
int main(int argc, char const *argv[]){
	// Initialize the MPI environment
  	MPI_Init(&argc, &argv);
  	//
  	int nprocs;
  	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  	//
  	int rank;
  	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  	// read AA from files named from proc1 to proc$n
	return 0;
}