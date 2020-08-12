
#if !defined(PETSCVIEWERHDF5_H)
#define PETSCVIEWERHDF5_H

#include <petscviewer.h>

#if defined(PETSC_HAVE_HDF5)
#include <hdf5.h>
#if !defined(H5_VERSION_GE)
/* H5_VERSION_GE was introduced in HDF5 1.8.7, we support >= 1.8.0 */
/* So beware this will automatically 0 for HDF5 1.8.0 - 1.8.6 */
#define H5_VERSION_GE(a,b,c) 0
#endif

PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetFileId(PetscViewer,hid_t*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5OpenGroup(PetscViewer, hid_t *, hid_t *);

/* On 32 bit systems HDF5 is limited by size of integer, because hsize_t is defined as size_t */
#define PETSC_HDF5_INT_MAX  (~(hsize_t)0)

/* As per https://portal.hdfgroup.org/display/HDF5/Chunking+in+HDF5, max. chunk size is 4GB */
#define PETSC_HDF5_MAX_CHUNKSIZE 2147483647

PETSC_STATIC_INLINE PetscErrorCode PetscHDF5IntCast(PetscInt a,hsize_t *b)
{
  PetscFunctionBegin;
  if (a < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Cannot convert negative size");
#if defined(PETSC_USE_64BIT_INDICES) && (H5_SIZEOF_HSIZE_T == 4)
  if (a > PETSC_HDF5_INT_MAX) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Array too long for HDF5");
#endif
  *b =  (hsize_t)(a);
  PetscFunctionReturn(0);
}
PETSC_EXTERN PetscErrorCode PetscDataTypeToHDF5DataType(PetscDataType,hid_t*);
PETSC_EXTERN PetscErrorCode PetscHDF5DataTypeToPetscDataType(hid_t,PetscDataType*);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5HasObject(PetscViewer,PetscObject,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5ReadAttribute(PetscViewer,const char[],const char[],PetscDataType,void*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5WriteAttribute(PetscViewer,const char[],const char[],PetscDataType,const void*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5HasAttribute(PetscViewer, const char[], const char[], PetscBool *);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5ReadObjectAttribute(PetscViewer,PetscObject,const char[],PetscDataType,void*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5WriteObjectAttribute(PetscViewer,PetscObject,const char[],PetscDataType,const void*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5HasObjectAttribute(PetscViewer,PetscObject,const char[],PetscBool*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5WriteSDS(PetscViewer,float *,int,int *,int);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5Open(MPI_Comm,const char[],PetscFileMode,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5PushGroup(PetscViewer,const char []);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5PopGroup(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetGroup(PetscViewer, const char *[]);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5HasGroup(PetscViewer,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5IncrementTimestep(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5SetTimestep(PetscViewer,PetscInt);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetTimestep(PetscViewer,PetscInt*);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5SetBaseDimension2(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetBaseDimension2(PetscViewer,PetscBool*);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5SetSPOutput(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetSPOutput(PetscViewer,PetscBool*);

PETSC_EXTERN PetscErrorCode PetscViewerHDF5SetCollective(PetscViewer,PetscBool);
PETSC_EXTERN PetscErrorCode PetscViewerHDF5GetCollective(PetscViewer,PetscBool*);
#endif  /* defined(PETSC_HAVE_HDF5) */
#endif
