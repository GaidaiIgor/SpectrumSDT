
#ifndef __VIEWERHDF5IMPL_H
#define __VIEWERHDF5IMPL_H

#if defined(H5_VERSION)
#  error "viewerhdf5impl.h must be included *before* any other HDF5 headers"
#else
#  define H5_USE_18_API
#endif
#include <petscviewerhdf5.h>

#if defined(PETSC_HAVE_HDF5)

#define PetscStackCallHDF5(func,args) do {                        \
    herr_t _status;                                               \
    PetscStackPush(#func);_status = func args;PetscStackPop; if (_status) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in HDF5 call %s() Status %d",#func,(int)_status); \
  } while (0)

#define PetscStackCallHDF5Return(ret,func,args) do {              \
    PetscStackPush(#func);ret = func args;PetscStackPop; if (ret < 0) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_LIB,"Error in HDF5 call %s() Status %d",#func,(int)ret); \
  } while (0)

typedef struct PetscViewerHDF5GroupList {
  const char       *name;
  struct PetscViewerHDF5GroupList *next;
} PetscViewerHDF5GroupList;

typedef struct {
  char          *filename;
  PetscFileMode btype;
  hid_t         file_id;
  hid_t         dxpl_id;   /* H5P_DATASET_XFER property list controlling raw data transfer (read/write). Properties are modified using H5Pset_dxpl_* functions. */
  PetscInt      timestep;
  PetscViewerHDF5GroupList *groups;
  PetscBool     basedimension2;  /* save vectors and DMDA vectors with a dimension of at least 2 even if the bs/dof is 1 */
  PetscBool     spoutput;  /* write data in single precision even if PETSc is compiled with double precision PetscReal */
} PetscViewer_HDF5;

#endif
#endif
