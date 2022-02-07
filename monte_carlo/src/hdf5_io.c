#include <hdf5.h>
#include "../include/global.h"
  /* The input for the hdf5 configuration */ 
  /* supposedly will read all the position and */ 
  /* triangles info from a single file */

void hdf5_io_read_config(double *Pos, int *cmlist,
        int *node_nbr, int2 *bond_nbr, int *triangles,
        char input_file[]){

    hid_t   file_id, group_id,dataset_id;  /* identifiers */
    herr_t  status;

  /* Open an existing file. */
  file_id = H5Fopen(input_file, H5F_ACC_RDONLY, H5P_DEFAULT);

  dataset_id = H5Dopen(file_id, "pos", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT,Pos);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "cumu_list", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, cmlist);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "node_nbr", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, node_nbr);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "bond_nbr", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, bond_nbr);
  status = H5Dclose(dataset_id);

  dataset_id = H5Dopen(file_id, "triangles", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, 
          H5S_ALL, H5S_ALL, H5P_DEFAULT, triangles);
  status = H5Dclose(dataset_id);

  status = H5Fclose(file_id);
}