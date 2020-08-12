#if !defined(__SFBASIC_H)
#define __SFBASIC_H

#include <petsc/private/sfimpl.h> /*I "petscsf.h" I*/

typedef struct _n_PetscSFLink* PetscSFLink;

#define SFBASICHEADER \
  PetscMPIInt      niranks;         /* Number of incoming ranks (ranks accessing my roots) */                                      \
  PetscMPIInt      ndiranks;        /* Number of incoming ranks (ranks accessing my roots) in distinguished set */                 \
  PetscMPIInt      *iranks;         /* Array of ranks that reference my roots */                                                   \
  PetscInt         itotal;          /* Total number of graph edges referencing my roots */                                         \
  PetscInt         *ioffset;        /* Array of length niranks+1 holding offset in irootloc[] for each rank */                     \
  PetscInt         *irootloc;       /* Incoming roots referenced by ranks starting at ioffset[rank] */                             \
  PetscInt         *irootloc_d[2];  /* A copy of irootloc[local/remote] in device memory if needed */                              \
  PetscInt         rootbuflen[2];   /* Length (in unit) of root buffers, in layout of [PETSCSF_LOCAL/REMOTE] */                    \
  PetscBool        rootcontig[2];   /* True means the local/remote segments of indices in irootloc[] are contiguous ... */         \
  PetscInt         rootstart[2];    /* ... and start from rootstart[0] and rootstart[1] respectively */                            \
  PetscSFPackOpt   rootpackopt[2];  /* Pack optimization plans based on patterns in irootloc[]. NULL for no optimizations */       \
  PetscSFPackOpt   rootpackopt_d[2];/* Copy of rootpackopt[] on device if needed */                                                \
  PetscBool        rootdups[2];     /* Indices of roots in irootloc[local/remote] have dups. Used for data-race test */            \
  PetscInt         nrootreqs;       /* Number of MPI reqests */                                                                    \
  PetscSFLink      avail;           /* One or more entries per MPI Datatype, lazily constructed */                                 \
  PetscSFLink      inuse            /* Buffers being used for transactions that have not yet completed */

typedef struct {
  SFBASICHEADER;
} PetscSF_Basic;

PETSC_STATIC_INLINE PetscErrorCode PetscSFGetRootInfo_Basic(PetscSF sf,PetscInt *nrootranks,PetscInt *ndrootranks,const PetscMPIInt **rootranks,const PetscInt **rootoffset,const PetscInt **rootloc)
{
  PetscSF_Basic *bas = (PetscSF_Basic*)sf->data;

  PetscFunctionBegin;
  if (nrootranks)  *nrootranks  = bas->niranks;
  if (ndrootranks) *ndrootranks = bas->ndiranks;
  if (rootranks)   *rootranks   = bas->iranks;
  if (rootoffset)  *rootoffset  = bas->ioffset;
  if (rootloc)     *rootloc     = bas->irootloc;
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode PetscSFGetLeafInfo_Basic(PetscSF sf,PetscInt *nleafranks,PetscInt *ndleafranks,const PetscMPIInt **leafranks,const PetscInt **leafoffset,const PetscInt **leafloc,const PetscInt **leafrremote)
{
  PetscFunctionBegin;
  if (nleafranks)  *nleafranks  = sf->nranks;
  if (ndleafranks) *ndleafranks = sf->ndranks;
  if (leafranks)   *leafranks   = sf->ranks;
  if (leafoffset)  *leafoffset  = sf->roffset;
  if (leafloc)     *leafloc     = sf->rmine;
  if (leafrremote) *leafrremote = sf->rremote;
  PetscFunctionReturn(0);
}

PETSC_INTERN PetscErrorCode PetscSFSetUp_Basic(PetscSF);
PETSC_INTERN PetscErrorCode PetscSFView_Basic(PetscSF,PetscViewer);
PETSC_INTERN PetscErrorCode PetscSFReset_Basic(PetscSF);
PETSC_INTERN PetscErrorCode PetscSFDestroy_Basic(PetscSF);
PETSC_INTERN PetscErrorCode PetscSFBcastAndOpEnd_Basic  (PetscSF,MPI_Datatype,const void*,void*,MPI_Op);
PETSC_INTERN PetscErrorCode PetscSFReduceEnd_Basic      (PetscSF,MPI_Datatype,const void*,void*,MPI_Op);
PETSC_INTERN PetscErrorCode PetscSFFetchAndOpBegin_Basic(PetscSF,MPI_Datatype,PetscMemType,void*,PetscMemType,const void*,void*,MPI_Op);
PETSC_INTERN PetscErrorCode PetscSFCreateEmbeddedSF_Basic(PetscSF,PetscInt,const PetscInt*,PetscSF*);
PETSC_INTERN PetscErrorCode PetscSFGetLeafRanks_Basic(PetscSF,PetscInt*,const PetscMPIInt**,const PetscInt**,const PetscInt**);
#endif
