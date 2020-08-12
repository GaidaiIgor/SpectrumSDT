
#include <petscfe.h>
#include <petscdmplex.h>
#include <petsc/private/hashmap.h>
#include <petsc/private/dmpleximpl.h>
#include <petsc/private/petscfeimpl.h>

const char help[] = "Test PETSCDUALSPACELAGRANGE\n";

typedef struct _PetscHashLagKey
{
  PetscInt  dim;
  PetscInt  order;
  PetscInt  formDegree;
  PetscBool trimmed;
  PetscInt  tensor;
  PetscBool continuous;
} PetscHashLagKey;

#define PetscHashLagKeyHash(key) \
  PetscHashCombine(PetscHashCombine(PetscHashCombine(PetscHashInt((key).dim), \
                                                     PetscHashInt((key).order)), \
                                    PetscHashInt((key).formDegree)), \
                   PetscHashCombine(PetscHashCombine(PetscHashInt((key).trimmed), \
                                                     PetscHashInt((key).tensor)), \
                                    PetscHashInt((key).continuous)))

#define PetscHashLagKeyEqual(k1,k2) \
  (((k1).dim == (k2).dim) ? \
   ((k1).order == (k2).order) ? \
   ((k1).formDegree == (k2).formDegree) ? \
   ((k1).trimmed == (k2).trimmed) ? \
   ((k1).tensor == (k2).tensor) ? \
   ((k1).continuous == (k2).continuous) : 0 : 0 : 0 : 0 : 0)

PETSC_HASH_MAP(HashLag, PetscHashLagKey, PetscInt, PetscHashLagKeyHash, PetscHashLagKeyEqual, 0)

static PetscErrorCode ExpectedNumDofs_Total(PetscInt dim, PetscInt order, PetscInt formDegree, PetscBool trimmed, PetscInt tensor, PetscInt nCopies, PetscInt *nDofs);
static PetscErrorCode ExpectedNumDofs_Interior(PetscInt dim, PetscInt order, PetscInt formDegree, PetscBool trimmed, PetscInt tensor, PetscInt nCopies, PetscInt *nDofs);

static PetscErrorCode ExpectedNumDofs_Total(PetscInt dim, PetscInt order, PetscInt formDegree, PetscBool trimmed, PetscInt tensor, PetscInt nCopies, PetscInt *nDofs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  formDegree = PetscAbsInt(formDegree);
  /* see femtable.org for the source of most of these values */
  *nDofs = -1;
  if (tensor == 0) { /* simplex spaces */
    if (trimmed) {
      PetscInt rnchooserk;
      PetscInt rkm1choosek;

      ierr = PetscDTBinomialInt(order + dim, order + formDegree, &rnchooserk);CHKERRQ(ierr);
      ierr = PetscDTBinomialInt(order + formDegree - 1, formDegree, &rkm1choosek);CHKERRQ(ierr);
      *nDofs = rnchooserk * rkm1choosek * nCopies;
    } else {
      PetscInt rnchooserk;
      PetscInt rkchoosek;

      ierr = PetscDTBinomialInt(order + dim, order + formDegree, &rnchooserk);CHKERRQ(ierr);
      ierr = PetscDTBinomialInt(order + formDegree, formDegree, &rkchoosek);CHKERRQ(ierr);
      *nDofs = rnchooserk * rkchoosek * nCopies;
    }
  } else if (tensor == 1) { /* hypercubes */
    if (trimmed) {
      PetscInt nchoosek;
      PetscInt rpowk, rp1pownmk;

      ierr = PetscDTBinomialInt(dim, formDegree, &nchoosek);CHKERRQ(ierr);
      rpowk = PetscPowInt(order, formDegree);CHKERRQ(ierr);
      rp1pownmk = PetscPowInt(order + 1, dim - formDegree);CHKERRQ(ierr);
      *nDofs = nchoosek * rpowk * rp1pownmk * nCopies;
    } else {
      PetscInt nchoosek;
      PetscInt rp1pown;

      ierr = PetscDTBinomialInt(dim, formDegree, &nchoosek);CHKERRQ(ierr);
      rp1pown = PetscPowInt(order + 1, dim);CHKERRQ(ierr);
      *nDofs = nchoosek * rp1pown * nCopies;
    }
  } else { /* prism */
    PetscInt tracek = 0;
    PetscInt tracekm1 = 0;
    PetscInt fiber0 = 0;
    PetscInt fiber1 = 0;

    if (formDegree < dim) {
      ierr = ExpectedNumDofs_Total(dim - 1, order, formDegree, trimmed, 0, 1, &tracek);CHKERRQ(ierr);
      ierr = ExpectedNumDofs_Total(1, order, 0, trimmed, 0, 1, &fiber0);CHKERRQ(ierr);
    }
    if (formDegree > 0) {
      ierr = ExpectedNumDofs_Total(dim - 1, order, formDegree - 1, trimmed, 0, 1, &tracekm1);CHKERRQ(ierr);
      ierr = ExpectedNumDofs_Total(1, order, 1, trimmed, 0, 1, &fiber1);CHKERRQ(ierr);
    }
    *nDofs = (tracek * fiber0 + tracekm1 * fiber1) * nCopies;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode ExpectedNumDofs_Interior(PetscInt dim, PetscInt order, PetscInt formDegree, PetscBool trimmed,
                                               PetscInt tensor, PetscInt nCopies, PetscInt *nDofs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  formDegree = PetscAbsInt(formDegree);
  /* see femtable.org for the source of most of these values */
  *nDofs = -1;
  if (tensor == 0) { /* simplex spaces */
    if (trimmed) {
      if (order + formDegree > dim) {
        PetscInt eorder = order + formDegree - dim - 1;
        PetscInt eformDegree = dim - formDegree;
        PetscInt rnchooserk;
        PetscInt rkchoosek;

        ierr = PetscDTBinomialInt(eorder + dim, eorder + eformDegree, &rnchooserk);CHKERRQ(ierr);
        ierr = PetscDTBinomialInt(eorder + eformDegree, eformDegree, &rkchoosek);CHKERRQ(ierr);
        *nDofs = rnchooserk * rkchoosek * nCopies;
      } else {
        *nDofs = 0;
      }

    } else {
      if (order + formDegree > dim) {
        PetscInt eorder = order + formDegree - dim;
        PetscInt eformDegree = dim - formDegree;
        PetscInt rnchooserk;
        PetscInt rkm1choosek;

        ierr = PetscDTBinomialInt(eorder + dim, eorder + eformDegree, &rnchooserk);CHKERRQ(ierr);
        ierr = PetscDTBinomialInt(eorder + eformDegree - 1, eformDegree, &rkm1choosek);CHKERRQ(ierr);
        *nDofs = rnchooserk * rkm1choosek * nCopies;
      } else {
        *nDofs = 0;
      }
    }
  } else if (tensor == 1) { /* hypercubes */
    if (dim < 2) {
      ierr = ExpectedNumDofs_Interior(dim, order, formDegree, trimmed, 0, nCopies, nDofs);CHKERRQ(ierr);
    } else {
      PetscInt tracek = 0;
      PetscInt tracekm1 = 0;
      PetscInt fiber0 = 0;
      PetscInt fiber1 = 0;

      if (formDegree < dim) {
        ierr = ExpectedNumDofs_Interior(dim - 1, order, formDegree, trimmed, dim > 2 ? 1 : 0, 1, &tracek);CHKERRQ(ierr);
        ierr = ExpectedNumDofs_Interior(1, order, 0, trimmed, 0, 1, &fiber0);CHKERRQ(ierr);
      }
      if (formDegree > 0) {
        ierr = ExpectedNumDofs_Interior(dim - 1, order, formDegree - 1, trimmed, dim > 2 ? 1 : 0, 1, &tracekm1);CHKERRQ(ierr);
        ierr = ExpectedNumDofs_Interior(1, order, 1, trimmed, 0, 1, &fiber1);CHKERRQ(ierr);
      }
      *nDofs = (tracek * fiber0 + tracekm1 * fiber1) * nCopies;
    }
  } else { /* prism */
    PetscInt tracek = 0;
    PetscInt tracekm1 = 0;
    PetscInt fiber0 = 0;
    PetscInt fiber1 = 0;

    if (formDegree < dim) {
      ierr = ExpectedNumDofs_Interior(dim - 1, order, formDegree, trimmed, 0, 1, &tracek);CHKERRQ(ierr);
      ierr = ExpectedNumDofs_Interior(1, order, 0, trimmed, 0, 1, &fiber0);CHKERRQ(ierr);
    }
    if (formDegree > 0) {
      ierr = ExpectedNumDofs_Interior(dim - 1, order, formDegree - 1, trimmed, 0, 1, &tracekm1);CHKERRQ(ierr);
      ierr = ExpectedNumDofs_Interior(1, order, 1, trimmed, 0, 1, &fiber1);CHKERRQ(ierr);
    }
    *nDofs = (tracek * fiber0 + tracekm1 * fiber1) * nCopies;
  }
  PetscFunctionReturn(0);
}

PetscErrorCode testLagrange(PetscHashLag lagTable, DM K, PetscInt dim, PetscInt order, PetscInt formDegree, PetscBool trimmed, PetscInt tensor, PetscBool continuous, PetscInt nCopies)
{
  PetscDualSpace  sp;
  MPI_Comm        comm = PETSC_COMM_SELF;
  PetscInt        Nk;
  PetscHashLagKey key;
  PetscHashIter   iter;
  PetscBool       missing;
  PetscInt        spdim, spintdim, exspdim, exspintdim;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscDTBinomialInt(dim, PetscAbsInt(formDegree), &Nk);CHKERRQ(ierr);
  ierr = PetscDualSpaceCreate(comm, &sp);CHKERRQ(ierr);
  ierr = PetscDualSpaceSetType(sp, PETSCDUALSPACELAGRANGE);CHKERRQ(ierr);
  ierr = PetscDualSpaceSetDM(sp, K);CHKERRQ(ierr);
  ierr = PetscDualSpaceSetOrder(sp, order);CHKERRQ(ierr);
  ierr = PetscDualSpaceSetFormDegree(sp, formDegree);CHKERRQ(ierr);
  ierr = PetscDualSpaceSetNumComponents(sp, nCopies * Nk);CHKERRQ(ierr);
  ierr = PetscDualSpaceLagrangeSetContinuity(sp, continuous);CHKERRQ(ierr);
  ierr = PetscDualSpaceLagrangeSetTensor(sp, (PetscBool) tensor);CHKERRQ(ierr);
  ierr = PetscDualSpaceLagrangeSetTrimmed(sp, trimmed);CHKERRQ(ierr);
  ierr = PetscInfo7(NULL, "Input: dim %D, order %D, trimmed %D, tensor %D, continuous %D, formDegree %D, nCopies %D\n", dim, order, (PetscInt) trimmed, tensor, (PetscInt) continuous, formDegree, nCopies);CHKERRQ(ierr);
  ierr = ExpectedNumDofs_Total(dim, order, formDegree, trimmed, tensor, nCopies, &exspdim);CHKERRQ(ierr);
  if (continuous && dim > 0 && order > 0) {
    ierr = ExpectedNumDofs_Interior(dim, order, formDegree, trimmed, tensor, nCopies, &exspintdim);CHKERRQ(ierr);
  } else {
    exspintdim = exspdim;
  }
  ierr = PetscDualSpaceSetUp(sp);CHKERRQ(ierr);
  ierr = PetscDualSpaceGetDimension(sp, &spdim);CHKERRQ(ierr);
  ierr = PetscDualSpaceGetInteriorDimension(sp, &spintdim);CHKERRQ(ierr);
  if (spdim != exspdim) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Expected dual space dimension %D, got %D\n", exspdim, spdim);CHKERRQ(ierr);
  if (spintdim != exspintdim) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Expected dual space interior dimension %D, got %D\n", exspintdim, spintdim);CHKERRQ(ierr);
  key.dim = dim;
  key.formDegree = formDegree;
  ierr = PetscDualSpaceGetOrder(sp, &key.order);CHKERRQ(ierr);
  ierr = PetscDualSpaceLagrangeGetContinuity(sp, &key.continuous);CHKERRQ(ierr);
  if (tensor == 2) {
    key.tensor = 2;
  } else {
    PetscBool bTensor;

    ierr = PetscDualSpaceLagrangeGetTensor(sp, &bTensor);CHKERRQ(ierr);
    key.tensor = bTensor;
  }
  ierr = PetscDualSpaceLagrangeGetTrimmed(sp, &key.trimmed);CHKERRQ(ierr);
  ierr = PetscInfo4(NULL, "After setup:  order %D, trimmed %D, tensor %D, continuous %D\n", key.order, (PetscInt) key.trimmed, key.tensor, (PetscInt) key.continuous);CHKERRQ(ierr);
  ierr = PetscHashLagPut(lagTable, key, &iter, &missing);CHKERRQ(ierr);
  if (missing) {
    DMPolytopeType type;

    ierr = DMPlexGetCellType(K, 0, &type);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "New space: %s, order %D, trimmed %D, tensor %D, continuous %D, form degree %D\n", DMPolytopeTypes[type], order, (PetscInt) trimmed, tensor, (PetscInt) continuous, formDegree);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
    {
      PetscQuadrature intNodes, allNodes;
      Mat intMat, allMat;
      MatInfo info;
      PetscInt i, j, nodeIdxDim, nodeVecDim, nNodes;
      const PetscInt *nodeIdx;
      const PetscReal *nodeVec;

      PetscDualSpace_Lag *lag = (PetscDualSpace_Lag *) sp->data;

      ierr = PetscLagNodeIndicesGetData_Internal(lag->allNodeIndices, &nodeIdxDim, &nodeVecDim, &nNodes, &nodeIdx, &nodeVec);CHKERRQ(ierr);
      if (nodeVecDim != Nk) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Incorrect nodeVecDim");CHKERRQ(ierr);
      if (nNodes != spdim) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Incorrect nNodes");CHKERRQ(ierr);

      ierr = PetscDualSpaceGetAllData(sp, &allNodes, &allMat);CHKERRQ(ierr);

      ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "All nodes:\n");CHKERRQ(ierr);
      ierr = PetscViewerASCIIPushTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
      ierr = PetscQuadratureView(allNodes, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPopTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "All node indices:\n");CHKERRQ(ierr);
      for (i = 0; i < spdim; i++) {
        ierr = PetscPrintf(PETSC_COMM_SELF, "(");CHKERRQ(ierr);
        for (j = 0; j < nodeIdxDim; j++) {
          ierr = PetscPrintf(PETSC_COMM_SELF, " %D,", nodeIdx[i * nodeIdxDim + j]);CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "): [");CHKERRQ(ierr);
        for (j = 0; j < nodeVecDim; j++) {
          ierr = PetscPrintf(PETSC_COMM_SELF, " %g,", (double) nodeVec[i * nodeVecDim + j]);CHKERRQ(ierr);
        }
        ierr = PetscPrintf(PETSC_COMM_SELF, "]\n");CHKERRQ(ierr);
      }

      ierr = MatGetInfo(allMat, MAT_LOCAL, &info);CHKERRQ(ierr);
      ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "All matrix: %D nonzeros\n", (PetscInt) info.nz_used);CHKERRQ(ierr);

      ierr = PetscDualSpaceGetInteriorData(sp, &intNodes, &intMat);CHKERRQ(ierr);
      if (intMat && intMat != allMat) {
        PetscInt intNodeIdxDim, intNodeVecDim, intNnodes;
        const PetscInt *intNodeIdx;
        const PetscReal *intNodeVec;
        PetscBool same;

        ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "Interior nodes:\n");CHKERRQ(ierr);
        ierr = PetscViewerASCIIPushTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
        ierr = PetscQuadratureView(intNodes, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPopTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);

        ierr = MatGetInfo(intMat, MAT_LOCAL, &info);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "Interior matrix: %D nonzeros\n", (PetscInt) info.nz_used);CHKERRQ(ierr);
        ierr = PetscLagNodeIndicesGetData_Internal(lag->intNodeIndices, &intNodeIdxDim, &intNodeVecDim, &intNnodes, &intNodeIdx, &intNodeVec);CHKERRQ(ierr);
        if (intNodeIdxDim != nodeIdxDim || intNodeVecDim != nodeVecDim) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Interior node indices not the same shale as all node indices");
        if (intNnodes != spintdim) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Incorrect interior nNodes");CHKERRQ(ierr);
        ierr = PetscArraycmp(intNodeIdx, nodeIdx, nodeIdxDim * intNnodes, &same);CHKERRQ(ierr);
        if (!same) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Interior node indices not the same as start of all node indices");
        ierr = PetscArraycmp(intNodeVec, nodeVec, nodeVecDim * intNnodes, &same);CHKERRQ(ierr);
        if (!same) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Interior node vectors not the same as start of all node vectors");
      } else if (intMat) {
        ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "Interior data is the same as all data\n");CHKERRQ(ierr);
        if (intNodes != allNodes) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Interior nodes should be the same as all nodes");
        if (lag->intNodeIndices != lag->allNodeIndices) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Interior node indices should be the same as all node indices");
      }
    }
    if (dim <= 2 && spintdim) {
      PetscInt coneSize, o;

      ierr = DMPlexGetConeSize(K, 0, &coneSize);CHKERRQ(ierr);
      for (o = -coneSize; o < coneSize; o++) {
        Mat symMat;

        ierr = PetscDualSpaceCreateInteriorSymmetryMatrix_Lagrange(sp, o, &symMat);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPrintf(PETSC_VIEWER_STDOUT_SELF, "Interior node symmetry matrix for orientation %D:\n", o);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPushTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
        ierr = MatView(symMat, PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
        ierr = PetscViewerASCIIPopTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
        ierr = MatDestroy(&symMat);CHKERRQ(ierr);
      }
    }
    ierr = PetscViewerASCIIPopTab(PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  }
  ierr = PetscDualSpaceDestroy(&sp);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCreateReferenceWedge(MPI_Comm comm, DM *refdm)
{
  PetscInt       dim = 3;
  DM             rdm;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMCreate(comm, &rdm);CHKERRQ(ierr);
  ierr = DMSetType(rdm, DMPLEX);CHKERRQ(ierr);
  ierr = DMSetDimension(rdm, dim);CHKERRQ(ierr);
  {
    PetscInt    numPoints[4]         = {6, 9, 5, 1};
    PetscInt    coneSize[21]         = {5,
                                        3, 3,
                                        4, 4, 4,
                                        2, 2, 2, 2, 2, 2, 2, 2, 2,
                                        0, 0, 0, 0, 0, 0};
    PetscInt    cones[41]            = {1, 2, 3, 4, 5,
                                        6, 7, 8,
                                        9, 10, 11,
                                        8, 12, 9, 13,
                                        7, 14, 10, 12,
                                        6, 13, 11, 14,
                                        15, 16,  16, 17,  17, 15,
                                        18, 19,  19, 20,  20, 18,
                                        17, 19,  18, 15,  16, 20};
    PetscInt    coneOrientations[41] = {0, 0, 0, 0, 0,
                                        0, 0, 0,
                                        0, 0, 0,
                                        -2,  0, -2,  0,
                                        -2,  0, -2, -2,
                                        -2, -2, -2, -2,
                                        0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0,
                                        0, 0, 0, 0, 0, 0};
    PetscScalar vertexCoords[18]    = {-1.0, -1.0, -1.0,
                                       -1.0,  1.0, -1.0,
                                        1.0, -1.0, -1.0,
                                       -1.0, -1.0,  1.0,
                                        1.0, -1.0,  1.0,
                                       -1.0,  1.0,  1.0};

    ierr = DMPlexCreateFromDAG(rdm, 3, numPoints, coneSize, cones, coneOrientations, vertexCoords);CHKERRQ(ierr);
  }
  *refdm = rdm;
  PetscFunctionReturn(0);
}

int main (int argc, char **argv)
{
  PetscInt        dim;
  PetscHashLag    lagTable;
  PetscInt        tensorCell;
  PetscInt        order, ordermin, ordermax;
  PetscBool       continuous;
  PetscBool       trimmed;
  DM              dm;
  PetscErrorCode  ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, help);if (ierr) return ierr;
  dim = 3;
  tensorCell = 0;
  continuous = PETSC_FALSE;
  trimmed = PETSC_FALSE;
  ierr = PetscOptionsBegin(PETSC_COMM_WORLD,"","Options for PETSCDUALSPACELAGRANGE test","none");CHKERRQ(ierr);
  ierr = PetscOptionsRangeInt("-dim", "The spatial dimension","ex1.c",dim,&dim,NULL,0,3);CHKERRQ(ierr);
  ierr = PetscOptionsRangeInt("-tensor", "(0) simplex (1) hypercube (2) wedge","ex1.c",tensorCell,&tensorCell,NULL,0,2);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-continuous", "Whether the dual space has continuity","ex1.c",continuous,&continuous,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-trimmed", "Whether the dual space matches a trimmed polynomial space","ex1.c",trimmed,&trimmed,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();
  ierr = PetscHashLagCreate(&lagTable);CHKERRQ(ierr);

  if (tensorCell < 2) {
    ierr = DMPlexCreateReferenceCell(PETSC_COMM_SELF, dim, (PetscBool) !tensorCell, &dm);CHKERRQ(ierr);
  } else {
    ierr = DMPlexCreateReferenceWedge(PETSC_COMM_SELF, &dm);CHKERRQ(ierr);
  }
  ordermin = trimmed ? 1 : 0;
  ordermax = tensorCell == 2 ? 4 : tensorCell == 1 ? 3 : dim + 2;
  for (order = ordermin; order <= ordermax; order++) {
    PetscInt formDegree;

    for (formDegree = PetscMin(0,-dim+1); formDegree <= dim; formDegree++) {
      PetscInt nCopies;

      for (nCopies = 1; nCopies <= 3; nCopies++) {
        ierr = testLagrange(lagTable, dm, dim, order, formDegree, trimmed, (PetscBool) tensorCell, continuous, nCopies);CHKERRQ(ierr);
      }
    }
  }
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscHashLagDestroy(&lagTable);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

 test:
   suffix: 0
   args: -dim 0

 test:
   suffix: 1_discontinuous_full
   args: -dim 1 -continuous 0 -trimmed 0

 test:
   suffix: 1_continuous_full
   args: -dim 1 -continuous 1 -trimmed 0

 test:
   suffix: 2_simplex_discontinuous_full
   args: -dim 2 -tensor 0 -continuous 0 -trimmed 0

 test:
   suffix: 2_simplex_continuous_full
   args: -dim 2 -tensor 0 -continuous 1 -trimmed 0

 test:
   suffix: 2_tensor_discontinuous_full
   args: -dim 2 -tensor 1 -continuous 0 -trimmed 0

 test:
   suffix: 2_tensor_continuous_full
   args: -dim 2 -tensor 1 -continuous 1 -trimmed 0

 test:
   suffix: 3_simplex_discontinuous_full
   args: -dim 3 -tensor 0 -continuous 0 -trimmed 0

 test:
   suffix: 3_simplex_continuous_full
   args: -dim 3 -tensor 0 -continuous 1 -trimmed 0

 test:
   suffix: 3_tensor_discontinuous_full
   args: -dim 3 -tensor 1 -continuous 0 -trimmed 0

 test:
   suffix: 3_tensor_continuous_full
   args: -dim 3 -tensor 1 -continuous 1 -trimmed 0

 test:
   suffix: 3_wedge_discontinuous_full
   args: -dim 3 -tensor 2 -continuous 0 -trimmed 0

 test:
   suffix: 3_wedge_continuous_full
   args: -dim 3 -tensor 2 -continuous 1 -trimmed 0

 test:
   suffix: 1_discontinuous_trimmed
   args: -dim 1 -continuous 0 -trimmed 1

 test:
   suffix: 1_continuous_trimmed
   args: -dim 1 -continuous 1 -trimmed 1

 test:
   suffix: 2_simplex_discontinuous_trimmed
   args: -dim 2 -tensor 0 -continuous 0 -trimmed 1

 test:
   suffix: 2_simplex_continuous_trimmed
   args: -dim 2 -tensor 0 -continuous 1 -trimmed 1

 test:
   suffix: 2_tensor_discontinuous_trimmed
   args: -dim 2 -tensor 1 -continuous 0 -trimmed 1

 test:
   suffix: 2_tensor_continuous_trimmed
   args: -dim 2 -tensor 1 -continuous 1 -trimmed 1

 test:
   suffix: 3_simplex_discontinuous_trimmed
   args: -dim 3 -tensor 0 -continuous 0 -trimmed 1

 test:
   suffix: 3_simplex_continuous_trimmed
   args: -dim 3 -tensor 0 -continuous 1 -trimmed 1

 test:
   suffix: 3_tensor_discontinuous_trimmed
   args: -dim 3 -tensor 1 -continuous 0 -trimmed 1

 test:
   suffix: 3_tensor_continuous_trimmed
   args: -dim 3 -tensor 1 -continuous 1 -trimmed 1

 test:
   suffix: 3_wedge_discontinuous_trimmed
   args: -dim 3 -tensor 2 -continuous 0 -trimmed 1

 test:
   suffix: 3_wedge_continuous_trimmed
   args: -dim 3 -tensor 2 -continuous 1 -trimmed 1

TEST*/
