static char help[] = "Poisson Problem in 2d and 3d with finite elements.\n\
We solve the Poisson problem in a rectangular\n\
domain, using a parallel unstructured mesh (DMPLEX) to discretize it.\n\
This example supports automatic convergence estimation\n\
and eventually adaptivity.\n\n\n";

#include <petscdmplex.h>
#include <petscsnes.h>
#include <petscds.h>
#include <petscconvest.h>

typedef struct {
  /* Domain and mesh definition */
  PetscInt  dim;               /* The topological mesh dimension */
  PetscBool simplex;           /* Simplicial mesh */
  PetscBool spectral;          /* Look at the spectrum along planes in the solution */
  PetscInt  cells[3];          /* The initial domain division */
  PetscBool shear;             /* Shear the domain */
  PetscBool adjoint;           /* Solve the adjoint problem */
} AppCtx;

static PetscErrorCode zero(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  *u = 0.0;
  return 0;
}

static PetscErrorCode trig_u(PetscInt dim, PetscReal time, const PetscReal x[], PetscInt Nc, PetscScalar *u, void *ctx)
{
  PetscInt d;
  *u = 0.0;
  for (d = 0; d < dim; ++d) *u += PetscSinReal(2.0*PETSC_PI*x[d]);
  return 0;
}

/* Compute integral of (residual of solution)*(adjoint solution - projection of adjoint solution) */
static void obj_error_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                        const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                        const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                        PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar obj[])
{
  obj[0] = a[aOff[0]]*(u[0] - a[aOff[1]]);
}

static void f0_trig_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                      const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                      const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                      PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) f0[0] += -4.0*PetscSqr(PETSC_PI)*PetscSinReal(2.0*PETSC_PI*x[d]);
}

static void f0_unity_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                       const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                       const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                       PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = 1.0;
}

static void f0_identityaux_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                             const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                             const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                             PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f0[])
{
  f0[0] = a[0];
}

static void f1_u(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                 const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                 const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                 PetscReal t, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar f1[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) f1[d] = u_x[d];
}

static void g3_uu(PetscInt dim, PetscInt Nf, PetscInt NfAux,
                  const PetscInt uOff[], const PetscInt uOff_x[], const PetscScalar u[], const PetscScalar u_t[], const PetscScalar u_x[],
                  const PetscInt aOff[], const PetscInt aOff_x[], const PetscScalar a[], const PetscScalar a_t[], const PetscScalar a_x[],
                  PetscReal t, PetscReal u_tShift, const PetscReal x[], PetscInt numConstants, const PetscScalar constants[], PetscScalar g3[])
{
  PetscInt d;
  for (d = 0; d < dim; ++d) g3[d*dim+d] = 1.0;
}

static PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscInt       n = 3;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  options->dim      = 2;
  options->cells[0] = 1;
  options->cells[1] = 1;
  options->cells[2] = 1;
  options->simplex  = PETSC_TRUE;
  options->shear    = PETSC_FALSE;
  options->spectral = PETSC_FALSE;
  options->adjoint  = PETSC_FALSE;

  ierr = PetscOptionsBegin(comm, "", "Poisson Problem Options", "DMPLEX");CHKERRQ(ierr);
  ierr = PetscOptionsInt("-dim", "The topological mesh dimension", "ex13.c", options->dim, &options->dim, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsIntArray("-cells", "The initial mesh division", "ex13.c", options->cells, &n, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-simplex", "Simplicial (true) or tensor (false) mesh", "ex13.c", options->simplex, &options->simplex, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-shear", "Shear the domain", "ex13.c", options->shear, &options->shear, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-spectral", "Look at the spectrum along planes of the solution", "ex13.c", options->spectral, &options->spectral, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsBool("-adjoint", "Solve the adjoint problem", "ex13.c", options->adjoint, &options->adjoint, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateSpectralPlanes(DM dm, PetscInt numPlanes, const PetscInt planeDir[], const PetscReal planeCoord[], AppCtx *user)
{
  PetscSection       coordSection;
  Vec                coordinates;
  const PetscScalar *coords;
  PetscInt           dim, p, vStart, vEnd, v;
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = DMGetCoordinateDim(dm, &dim);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
  ierr = DMGetCoordinateSection(dm, &coordSection);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coordinates, &coords);CHKERRQ(ierr);
  for (p = 0; p < numPlanes; ++p) {
    DMLabel label;
    char    name[PETSC_MAX_PATH_LEN];

    ierr = PetscSNPrintf(name, PETSC_MAX_PATH_LEN, "spectral_plane_%D", p);CHKERRQ(ierr);
    ierr = DMCreateLabel(dm, name);CHKERRQ(ierr);
    ierr = DMGetLabel(dm, name, &label);CHKERRQ(ierr);
    ierr = DMLabelAddStratum(label, 1);CHKERRQ(ierr);
    for (v = vStart; v < vEnd; ++v) {
      PetscInt off;

      ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
      if (PetscAbsReal(planeCoord[p] - PetscRealPart(coords[off+planeDir[p]])) < PETSC_SMALL) {
	ierr = DMLabelSetValue(label, v, 1);CHKERRQ(ierr);
      }
    }
  }
  ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode CreateMesh(MPI_Comm comm, AppCtx *user, DM *dm)
{
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* Create box mesh */
  ierr = DMPlexCreateBoxMesh(comm, user->dim, user->simplex, user->cells, NULL, NULL, NULL, PETSC_TRUE, dm);CHKERRQ(ierr);
  /* Distribute mesh over processes */
  {
    DM               dmDist = NULL;
    PetscPartitioner part;

    ierr = DMPlexGetPartitioner(*dm, &part);CHKERRQ(ierr);
    ierr = PetscPartitionerSetFromOptions(part);CHKERRQ(ierr);
    ierr = DMPlexDistribute(*dm, 0, NULL, &dmDist);CHKERRQ(ierr);
    if (dmDist) {
      ierr = DMDestroy(dm);CHKERRQ(ierr);
      *dm  = dmDist;
    }
  }
  /* TODO: This should be pulled into the library */
  {
    char      convType[256];
    PetscBool flg;

    ierr = PetscOptionsBegin(comm, "", "Mesh conversion options", "DMPLEX");CHKERRQ(ierr);
    ierr = PetscOptionsFList("-dm_plex_convert_type","Convert DMPlex to another format","ex12",DMList,DMPLEX,convType,256,&flg);CHKERRQ(ierr);
    ierr = PetscOptionsEnd();
    if (flg) {
      DM dmConv;

      ierr = DMConvert(*dm,convType,&dmConv);CHKERRQ(ierr);
      if (dmConv) {
        ierr = DMDestroy(dm);CHKERRQ(ierr);
        *dm  = dmConv;
      }
    }
  }
  if (user->shear) {ierr = DMPlexShearGeometry(*dm, DM_X, NULL);CHKERRQ(ierr);}
  /* TODO: This should be pulled into the library */
  ierr = DMLocalizeCoordinates(*dm);CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) *dm, "Mesh");CHKERRQ(ierr);
  ierr = DMSetApplicationContext(*dm, user);CHKERRQ(ierr);
  ierr = DMSetFromOptions(*dm);CHKERRQ(ierr);
  ierr = DMViewFromOptions(*dm, NULL, "-dm_view");CHKERRQ(ierr);
  /* TODO: Add a hierachical viewer */
  if (user->spectral) {
    PetscInt  planeDir[2]   = {0,  1};
    PetscReal planeCoord[2] = {0., 1.};

    ierr = CreateSpectralPlanes(*dm, 2, planeDir, planeCoord, user);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupPrimalProblem(DM dm, AppCtx *user)
{
  PetscDS        prob;
  const PetscInt id = 1;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);
  ierr = PetscDSSetResidual(prob, 0, f0_trig_u, f1_u);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, 0, 0, NULL, NULL, NULL, g3_uu);CHKERRQ(ierr);
  ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "wall", "marker", 0, 0, NULL, (void (*)(void)) trig_u, 1, &id, user);CHKERRQ(ierr);
  ierr = PetscDSSetExactSolution(prob, 0, trig_u, user);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupAdjointProblem(DM dm, AppCtx *user)
{
  PetscDS        prob;
  const PetscInt id = 1;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);
  ierr = PetscDSSetResidual(prob, 0, f0_unity_u, f1_u);CHKERRQ(ierr);
  ierr = PetscDSSetJacobian(prob, 0, 0, NULL, NULL, NULL, g3_uu);CHKERRQ(ierr);
  ierr = PetscDSAddBoundary(prob, DM_BC_ESSENTIAL, "wall", "marker", 0, 0, NULL, (void (*)(void)) zero, 1, &id, user);CHKERRQ(ierr);
  ierr = PetscDSSetObjective(prob, 0, obj_error_u);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupErrorProblem(DM dm, AppCtx *user)
{
  PetscDS        prob;
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  ierr = DMGetDS(dm, &prob);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode SetupDiscretization(DM dm, const char name[], PetscErrorCode (*setup)(DM, AppCtx *), AppCtx *user)
{
  DM             cdm = dm;
  PetscFE        fe;
  char           prefix[PETSC_MAX_PATH_LEN];
  PetscErrorCode ierr;

  PetscFunctionBeginUser;
  /* Create finite element */
  ierr = PetscSNPrintf(prefix, PETSC_MAX_PATH_LEN, "%s_", name);CHKERRQ(ierr);
  ierr = PetscFECreateDefault(PetscObjectComm((PetscObject) dm), user->dim, 1, user->simplex, name ? prefix : NULL, -1, &fe);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) fe, name);CHKERRQ(ierr);
  /* Set discretization and boundary conditions for each mesh */
  ierr = DMSetField(dm, 0, NULL, (PetscObject) fe);CHKERRQ(ierr);
  ierr = DMCreateDS(dm);CHKERRQ(ierr);
  ierr = (*setup)(dm, user);CHKERRQ(ierr);
  while (cdm) {
    ierr = DMCopyDisc(dm,cdm);CHKERRQ(ierr);
    /* TODO: Check whether the boundary of coarse meshes is marked */
    ierr = DMGetCoarseDM(cdm, &cdm);CHKERRQ(ierr);
  }
  ierr = PetscFEDestroy(&fe);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode ComputeSpectral(DM dm, Vec u, PetscInt numPlanes, const PetscInt planeDir[], const PetscReal planeCoord[], AppCtx *user)
{
  MPI_Comm           comm;
  PetscSection       coordSection, section;
  Vec                coordinates, uloc;
  const PetscScalar *coords, *array;
  PetscInt           p;
  PetscMPIInt        size, rank;
  PetscErrorCode     ierr;

  PetscFunctionBeginUser;
  ierr = PetscObjectGetComm((PetscObject) dm, &comm);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = DMGetLocalVector(dm, &uloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(dm, u, INSERT_VALUES, uloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(dm, u, INSERT_VALUES, uloc);CHKERRQ(ierr);
  ierr = DMPlexInsertBoundaryValues(dm, PETSC_TRUE, uloc, 0.0, NULL, NULL, NULL);CHKERRQ(ierr);
  ierr = VecViewFromOptions(uloc, NULL, "-sol_view");CHKERRQ(ierr);
  ierr = DMGetLocalSection(dm, &section);CHKERRQ(ierr);
  ierr = VecGetArrayRead(uloc, &array);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
  ierr = DMGetCoordinateSection(dm, &coordSection);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coordinates, &coords);CHKERRQ(ierr);
  for (p = 0; p < numPlanes; ++p) {
    DMLabel         label;
    char            name[PETSC_MAX_PATH_LEN];
    Mat             F;
    Vec             x, y;
    IS              stratum;
    PetscReal      *ray, *gray;
    PetscScalar    *rvals, *svals, *gsvals;
    PetscInt       *perm, *nperm;
    PetscInt        n, N, i, j, off, offu;
    const PetscInt *points;

    ierr = PetscSNPrintf(name, PETSC_MAX_PATH_LEN, "spectral_plane_%D", p);CHKERRQ(ierr);
    ierr = DMGetLabel(dm, name, &label);CHKERRQ(ierr);
    ierr = DMLabelGetStratumIS(label, 1, &stratum);CHKERRQ(ierr);
    ierr = ISGetLocalSize(stratum, &n);CHKERRQ(ierr);
    ierr = ISGetIndices(stratum, &points);CHKERRQ(ierr);
    ierr = PetscMalloc2(n, &ray, n, &svals);CHKERRQ(ierr);
    for (i = 0; i < n; ++i) {
      ierr = PetscSectionGetOffset(coordSection, points[i], &off);CHKERRQ(ierr);
      ierr = PetscSectionGetOffset(section, points[i], &offu);CHKERRQ(ierr);
      ray[i]   = PetscRealPart(coords[off+((planeDir[p]+1)%2)]);
      svals[i] = array[offu];
    }
    /* Gather the ray data to proc 0 */
    if (size > 1) {
      PetscMPIInt *cnt, *displs, p;

      ierr = PetscCalloc2(size, &cnt, size, &displs);CHKERRQ(ierr);
      ierr = MPI_Gather(&n, 1, MPIU_INT, cnt, 1, MPIU_INT, 0, comm);CHKERRQ(ierr);
      for (p = 1; p < size; ++p) displs[p] = displs[p-1] + cnt[p-1];
      N = displs[size-1] + cnt[size-1];
      ierr = PetscMalloc2(N, &gray, N, &gsvals);CHKERRQ(ierr);
      ierr = MPI_Gatherv(ray, n, MPIU_REAL, gray, cnt, displs, MPIU_REAL, 0, comm);CHKERRQ(ierr);
      ierr = MPI_Gatherv(svals, n, MPIU_SCALAR, gsvals, cnt, displs, MPIU_SCALAR, 0, comm);CHKERRQ(ierr);
      ierr = PetscFree2(cnt, displs);CHKERRQ(ierr);
    } else {
      N      = n;
      gray   = ray;
      gsvals = svals;
    }
    if (!rank) {
      /* Sort point along ray */
      ierr = PetscMalloc2(N, &perm, N, &nperm);CHKERRQ(ierr);
      for (i = 0; i < N; ++i) {perm[i] = i;}
      ierr = PetscSortRealWithPermutation(N, gray, perm);CHKERRQ(ierr);
      /* Count duplicates and squish mapping */
      nperm[0] = perm[0];
      for (i = 1, j = 1; i < N; ++i) {
        if (PetscAbsReal(gray[perm[i]] - gray[perm[i-1]]) > PETSC_SMALL) nperm[j++] = perm[i];
      }
      /* Create FFT structs */
      ierr = MatCreateFFT(PETSC_COMM_SELF, 1, &j, MATFFTW, &F);CHKERRQ(ierr);
      ierr = MatCreateVecs(F, &x, &y);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject) y, name);CHKERRQ(ierr);
      ierr = VecGetArray(x, &rvals);CHKERRQ(ierr);
      for (i = 0, j = 0; i < N; ++i) {
        if (i > 0 && PetscAbsReal(gray[perm[i]] - gray[perm[i-1]]) < PETSC_SMALL) continue;
        rvals[j] = gsvals[nperm[j]];
        ++j;
      }
      ierr = PetscFree2(perm, nperm);CHKERRQ(ierr);
      if (size > 1) {ierr = PetscFree2(gray, gsvals);CHKERRQ(ierr);}
      ierr = VecRestoreArray(x, &rvals);CHKERRQ(ierr);
      /* Do FFT along the ray */
      ierr = MatMult(F, x, y);CHKERRQ(ierr);
      /* Chop FFT */
      ierr = VecChop(y, PETSC_SMALL);CHKERRQ(ierr);
      ierr = VecViewFromOptions(x, NULL, "-real_view");CHKERRQ(ierr);
      ierr = VecViewFromOptions(y, NULL, "-fft_view");CHKERRQ(ierr);
      ierr = VecDestroy(&x);CHKERRQ(ierr);
      ierr = VecDestroy(&y);CHKERRQ(ierr);
      ierr = MatDestroy(&F);CHKERRQ(ierr);
    }
    ierr = ISRestoreIndices(stratum, &points);CHKERRQ(ierr);
    ierr = ISDestroy(&stratum);CHKERRQ(ierr);
    ierr = PetscFree2(ray, svals);CHKERRQ(ierr);
  }
  ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(uloc, &array);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(dm, &uloc);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

int main(int argc, char **argv)
{
  DM             dm;   /* Problem specification */
  SNES           snes; /* Nonlinear solver */
  Vec            u;    /* Solutions */
  AppCtx         user; /* User-defined work context */
  PetscErrorCode ierr;

  ierr = PetscInitialize(&argc, &argv, NULL,help);if (ierr) return ierr;
  ierr = ProcessOptions(PETSC_COMM_WORLD, &user);CHKERRQ(ierr);
  /* Primal system */
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);CHKERRQ(ierr);
  ierr = CreateMesh(PETSC_COMM_WORLD, &user, &dm);CHKERRQ(ierr);
  ierr = SNESSetDM(snes, dm);CHKERRQ(ierr);
  ierr = SetupDiscretization(dm, "potential", SetupPrimalProblem, &user);CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(dm, &u);CHKERRQ(ierr);
  ierr = VecSet(u, 0.0);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) u, "potential");CHKERRQ(ierr);
  ierr = DMPlexSetSNESLocalFEM(dm, &user, &user, &user);CHKERRQ(ierr);
  ierr = SNESSetFromOptions(snes);CHKERRQ(ierr);
  ierr = SNESSolve(snes, NULL, u);CHKERRQ(ierr);
  ierr = SNESGetSolution(snes, &u);CHKERRQ(ierr);
  ierr = VecViewFromOptions(u, NULL, "-potential_view");CHKERRQ(ierr);
  if (user.spectral) {
    PetscInt  planeDir[2]   = {0,  1};
    PetscReal planeCoord[2] = {0., 1.};

    ierr = ComputeSpectral(dm, u, 2, planeDir, planeCoord, &user);CHKERRQ(ierr);
  }
  /* Adjoint system */
  if (user.adjoint) {
    DM   dmAdj;
    SNES snesAdj;
    Vec  uAdj;

    ierr = SNESCreate(PETSC_COMM_WORLD, &snesAdj);CHKERRQ(ierr);
    ierr = PetscObjectSetOptionsPrefix((PetscObject) snesAdj, "adjoint_");CHKERRQ(ierr);
    ierr = DMClone(dm, &dmAdj);CHKERRQ(ierr);
    ierr = SNESSetDM(snesAdj, dmAdj);CHKERRQ(ierr);
    ierr = SetupDiscretization(dmAdj, "adjoint", SetupAdjointProblem, &user);CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(dmAdj, &uAdj);CHKERRQ(ierr);
    ierr = VecSet(uAdj, 0.0);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) uAdj, "adjoint");CHKERRQ(ierr);
    ierr = DMPlexSetSNESLocalFEM(dmAdj, &user, &user, &user);CHKERRQ(ierr);
    ierr = SNESSetFromOptions(snesAdj);CHKERRQ(ierr);
    ierr = SNESSolve(snesAdj, NULL, uAdj);CHKERRQ(ierr);
    ierr = SNESGetSolution(snesAdj, &uAdj);CHKERRQ(ierr);
    ierr = VecViewFromOptions(uAdj, NULL, "-adjoint_view");CHKERRQ(ierr);
    /* Error representation */
    {
      DM        dmErr, dmErrAux, dms[2];
      Vec       errorEst, errorL2, uErr, uErrLoc, uAdjLoc, uAdjProj;
      IS       *subis;
      PetscReal errorEstTot, errorL2Norm, errorL2Tot;
      PetscInt  N, i;
      PetscErrorCode (*funcs[1])(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *) = {trig_u};
      void (*identity[1])(PetscInt, PetscInt, PetscInt,
                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                          const PetscInt[], const PetscInt[], const PetscScalar[], const PetscScalar[], const PetscScalar[],
                          PetscReal, const PetscReal[], PetscInt, const PetscScalar[], PetscScalar[]) = {f0_identityaux_u};
      void            *ctxs[1] = {0};

      ctxs[0] = &user;
      ierr = DMClone(dm, &dmErr);CHKERRQ(ierr);
      ierr = SetupDiscretization(dmErr, "error", SetupErrorProblem, &user);CHKERRQ(ierr);
      ierr = DMGetGlobalVector(dmErr, &errorEst);CHKERRQ(ierr);
      ierr = DMGetGlobalVector(dmErr, &errorL2);CHKERRQ(ierr);
      /*   Compute auxiliary data (solution and projection of adjoint solution) */
      ierr = DMGetLocalVector(dmAdj, &uAdjLoc);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(dmAdj, uAdj, INSERT_VALUES, uAdjLoc);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dmAdj, uAdj, INSERT_VALUES, uAdjLoc);CHKERRQ(ierr);
      ierr = DMGetGlobalVector(dm, &uAdjProj);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dm, "dmAux", (PetscObject) dmAdj);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dm, "A", (PetscObject) uAdjLoc);CHKERRQ(ierr);
      ierr = DMProjectField(dm, 0.0, u, identity, INSERT_VALUES, uAdjProj);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dm, "dmAux", NULL);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dm, "A", NULL);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dmAdj, &uAdjLoc);CHKERRQ(ierr);
      /*   Attach auxiliary data */
      dms[0] = dm; dms[1] = dm;
      ierr = DMCreateSuperDM(dms, 2, &subis, &dmErrAux);CHKERRQ(ierr);
      if (0) {
        PetscSection sec;

        ierr = DMGetLocalSection(dms[0], &sec);CHKERRQ(ierr);
        ierr = PetscSectionView(sec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        ierr = DMGetLocalSection(dms[1], &sec);CHKERRQ(ierr);
        ierr = PetscSectionView(sec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
        ierr = DMGetLocalSection(dmErrAux, &sec);CHKERRQ(ierr);
        ierr = PetscSectionView(sec, PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
      }
      ierr = DMViewFromOptions(dmErrAux, NULL, "-dm_err_view");CHKERRQ(ierr);
      ierr = ISViewFromOptions(subis[0], NULL, "-super_is_view");CHKERRQ(ierr);
      ierr = ISViewFromOptions(subis[1], NULL, "-super_is_view");CHKERRQ(ierr);
      ierr = DMGetGlobalVector(dmErrAux, &uErr);CHKERRQ(ierr);
      ierr = VecViewFromOptions(u, NULL, "-map_vec_view");CHKERRQ(ierr);
      ierr = VecViewFromOptions(uAdjProj, NULL, "-map_vec_view");CHKERRQ(ierr);
      ierr = VecViewFromOptions(uErr, NULL, "-map_vec_view");CHKERRQ(ierr);
      ierr = VecISCopy(uErr, subis[0], SCATTER_FORWARD, u);CHKERRQ(ierr);
      ierr = VecISCopy(uErr, subis[1], SCATTER_FORWARD, uAdjProj);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(dm, &uAdjProj);CHKERRQ(ierr);
      for (i = 0; i < 2; ++i) {ierr = ISDestroy(&subis[i]);CHKERRQ(ierr);}
      ierr = PetscFree(subis);CHKERRQ(ierr);
      ierr = DMGetLocalVector(dmErrAux, &uErrLoc);CHKERRQ(ierr);
      ierr = DMGlobalToLocalBegin(dm, uErr, INSERT_VALUES, uErrLoc);CHKERRQ(ierr);
      ierr = DMGlobalToLocalEnd(dm, uErr, INSERT_VALUES, uErrLoc);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(dmErrAux, &uErr);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dmAdj, "dmAux", (PetscObject) dmErrAux);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dmAdj, "A", (PetscObject) uErrLoc);CHKERRQ(ierr);
      /*   Compute cellwise error estimate */
      ierr = VecSet(errorEst, 0.0);CHKERRQ(ierr);
      ierr = DMPlexComputeCellwiseIntegralFEM(dmAdj, uAdj, errorEst, &user);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dmAdj, "dmAux", NULL);CHKERRQ(ierr);
      ierr = PetscObjectCompose((PetscObject) dmAdj, "A", NULL);CHKERRQ(ierr);
      ierr = DMRestoreLocalVector(dmErrAux, &uErrLoc);CHKERRQ(ierr);
      ierr = DMDestroy(&dmErrAux);CHKERRQ(ierr);
      /*   Plot cellwise error vector */
      ierr = VecViewFromOptions(errorEst, NULL, "-error_view");CHKERRQ(ierr);
      /*   Compute ratio of estimate (sum over cells) with actual L_2 error */
      ierr = DMComputeL2Diff(dm, 0.0, funcs, ctxs, u, &errorL2Norm);CHKERRQ(ierr);
      ierr = DMPlexComputeL2DiffVec(dm, 0.0, funcs, ctxs, u, errorL2);CHKERRQ(ierr);
      ierr = VecViewFromOptions(errorL2, NULL, "-l2_error_view");CHKERRQ(ierr);
      ierr = VecNorm(errorL2,  NORM_INFINITY, &errorL2Tot);CHKERRQ(ierr);
      ierr = VecNorm(errorEst, NORM_INFINITY, &errorEstTot);CHKERRQ(ierr);
      ierr = VecGetSize(errorEst, &N);CHKERRQ(ierr);
      ierr = VecPointwiseDivide(errorEst, errorEst, errorL2);CHKERRQ(ierr);
      ierr = PetscObjectSetName((PetscObject) errorEst, "Error ratio");CHKERRQ(ierr);
      ierr = VecViewFromOptions(errorEst, NULL, "-error_ratio_view");CHKERRQ(ierr);
      ierr = PetscPrintf(PETSC_COMM_WORLD, "N: %D L2 error: %g Error Ratio: %g/%g = %g\n", N, (double) errorL2Norm, (double) errorEstTot, (double) PetscSqrtReal(errorL2Tot), (double) errorEstTot/PetscSqrtReal(errorL2Tot));CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(dmErr, &errorEst);CHKERRQ(ierr);
      ierr = DMRestoreGlobalVector(dmErr, &errorL2);CHKERRQ(ierr);
      ierr = DMDestroy(&dmErr);CHKERRQ(ierr);
    }
    ierr = DMDestroy(&dmAdj);CHKERRQ(ierr);
    ierr = VecDestroy(&uAdj);CHKERRQ(ierr);
    ierr = SNESDestroy(&snesAdj);CHKERRQ(ierr);
  }
  /* Cleanup */
  ierr = VecDestroy(&u);CHKERRQ(ierr);
  ierr = SNESDestroy(&snes);CHKERRQ(ierr);
  ierr = DMDestroy(&dm);CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

  test:
    suffix: 2d_p1_0
    requires: triangle
    args: -potential_petscspace_degree 1 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_p1_scalable
    requires: triangle long_runtime
    args: -potential_petscspace_order 1 -dm_refine 3 -num_refine 3 -snes_convergence_estimate \
      -ksp_type cg -ksp_rtol 1.e-11 -ksp_norm_type unpreconditioned \
      -pc_type gamg \
        -pc_gamg_type agg -pc_gamg_agg_nsmooths 1 \
        -pc_gamg_coarse_eq_limit 1000 \
        -pc_gamg_square_graph 1 \
        -pc_gamg_threshold 0.05 \
        -pc_gamg_threshold_scale .0 \
      -mg_levels_ksp_type chebyshev \
        -mg_levels_ksp_max_it 1 \
        -mg_levels_esteig_ksp_type cg \
        -mg_levels_esteig_ksp_max_it 10 \
        -mg_levels_ksp_chebyshev_esteig 0,0.05,0,1.05 \
        -mg_levels_pc_type jacobi \
      -matptap_via scalable
  test:
    suffix: 2d_p1_gmg_vcycle
    requires: triangle
    args: -potential_petscspace_degree 1 -cells 2,2 -dm_refine_hierarchy 2 -convest_num_refine 2 -snes_convergence_estimate \
          -ksp_rtol 5e-10 -pc_type mg \
            -mg_levels_ksp_max_it 1 \
            -mg_levels_esteig_ksp_type cg \
            -mg_levels_esteig_ksp_max_it 10 \
            -mg_levels_ksp_chebyshev_esteig 0,0.05,0,1.05 \
            -mg_levels_pc_type jacobi
  test:
    suffix: 2d_p1_gmg_fcycle
    requires: triangle
    args: -potential_petscspace_degree 1 -cells 2,2 -dm_refine_hierarchy 2 -convest_num_refine 2 -snes_convergence_estimate \
          -ksp_rtol 5e-10 -pc_type mg -pc_mg_type full \
            -mg_levels_ksp_max_it 2 \
            -mg_levels_esteig_ksp_type cg \
            -mg_levels_esteig_ksp_max_it 10 \
            -mg_levels_ksp_chebyshev_esteig 0,0.05,0,1.05 \
            -mg_levels_pc_type jacobi
  test:
    suffix: 2d_p2_0
    requires: triangle
    args: -potential_petscspace_degree 2 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_p3_0
    requires: triangle
    args: -potential_petscspace_degree 3 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_q1_0
    args: -simplex 0 -potential_petscspace_degree 1 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_q1_1
    args: -simplex 0 -shear -potential_petscspace_degree 1 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_q2_0
    args: -simplex 0 -potential_petscspace_degree 2 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_q2_1
    args: -simplex 0 -shear -potential_petscspace_degree 2 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_q3_0
    args: -simplex 0 -potential_petscspace_degree 3 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 2d_q3_1
    args: -simplex 0 -shear -potential_petscspace_degree 3 -dm_refine 2 -convest_num_refine 3 -snes_convergence_estimate
  test:
    suffix: 3d_p1_0
    requires: ctetgen
    args: -dim 3 -cells 2,2,2 -potential_petscspace_degree 1 -convest_num_refine 2 -snes_convergence_estimate
  test:
    suffix: 3d_p2_0
    requires: ctetgen
    args: -dim 3 -cells 2,2,2 -potential_petscspace_degree 2 -convest_num_refine 2 -snes_convergence_estimate
  test:
    suffix: 3d_p3_0
    requires: ctetgen
    timeoutfactor: 2
    args: -dim 3 -cells 2,2,2 -potential_petscspace_degree 3 -convest_num_refine 2 -snes_convergence_estimate
  test:
    suffix: 3d_q1_0
    args: -dim 3 -simplex 0 -potential_petscspace_degree 1 -dm_refine 1 -convest_num_refine 2 -snes_convergence_estimate
  test:
    suffix: 3d_q2_0
    args: -dim 3 -simplex 0 -potential_petscspace_degree 2 -dm_refine 1 -convest_num_refine 2 -snes_convergence_estimate
  test:
    suffix: 3d_q3_0
    args: -dim 3 -simplex 0 -potential_petscspace_degree 3 -convest_num_refine 2 -snes_convergence_estimate
  test:
    suffix: 2d_p1_spectral_0
    requires: triangle fftw !complex
    args: -potential_petscspace_degree 1 -dm_refine 6 -spectral -fft_view
  test:
    suffix: 2d_p1_spectral_1
    requires: triangle fftw !complex
    nsize: 2
    args: -potential_petscspace_degree 1 -dm_refine 2 -spectral -fft_view
  test:
    suffix: 2d_p1_adj_0
    requires: triangle
    args: -potential_petscspace_degree 1 -dm_refine 2 -adjoint -adjoint_petscspace_degree 1 -error_petscspace_degree 0

TEST*/
