#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include <petsc/private/petscfeimpl.h>  /* For PetscFEInterpolate_Static() */
#include <petscsf.h>

const char * const DMPlexCellRefinerTypes[] = {"Regular", "ToBox", "ToSimplex", "DMPlexCellRefinerTypes", "DM_REFINER_", 0};

/*
  Note that j and invj are non-square:
         v0 + j x_face = x_cell
    invj (x_cell - v0) = x_face
*/
static PetscErrorCode DMPlexCellRefinerGetAffineFaceTransforms_Regular(DMPlexCellRefiner cr, DMPolytopeType ct, PetscInt *Nf, PetscReal *v0[], PetscReal *J[], PetscReal *invJ[], PetscReal *detJ[])
{
  /*
   2
   |\
   | \
   |  \
   |   \
   |    \
   |     \
   |      \
   2       1
   |        \
   |         \
   |          \
   0---0-------1
   v0[Nf][dc]:       3 x 2
   J[Nf][df][dc]:    3 x 1 x 2
   invJ[Nf][dc][df]: 3 x 2 x 1
   detJ[Nf]:         3
   */
  static PetscReal tri_v0[]   = {0.0, -1.0,  0.0, 0.0,  -1.0,  0.0};
  static PetscReal tri_J[]    = {1.0, 0.0,  -1.0, 1.0,   0.0, -1.0};
  static PetscReal tri_invJ[] = {1.0, 0.0,  -0.5, 0.5,   0.0, -1.0};
  static PetscReal tri_detJ[] = {1.0,  1.414213562373095,  1.0};
  /*
   3---------2---------2
   |                   |
   |                   |
   |                   |
   3                   1
   |                   |
   |                   |
   |                   |
   0---------0---------1

   v0[Nf][dc]:       4 x 2
   J[Nf][df][dc]:    4 x 1 x 2
   invJ[Nf][dc][df]: 4 x 2 x 1
   detJ[Nf]:         4
   */
  static PetscReal quad_v0[]   = {0.0, -1.0,  1.0, 0.0,   0.0, 1.0  -1.0,  0.0};
  static PetscReal quad_J[]    = {1.0, 0.0,   0.0, 1.0,  -1.0, 0.0,  0.0, -1.0};
  static PetscReal quad_invJ[] = {1.0, 0.0,   0.0, 1.0,  -1.0, 0.0,  0.0, -1.0};
  static PetscReal quad_detJ[] = {1.0,  1.0,  1.0,  1.0};

  PetscFunctionBegin;
  switch (ct) {
  case DM_POLYTOPE_TRIANGLE:      if (Nf) *Nf = 3; if (v0) *v0 = tri_v0;  if (J) *J = tri_J;  if (invJ) *invJ = tri_invJ;  if (detJ) *detJ = tri_detJ;  break;
  case DM_POLYTOPE_QUADRILATERAL: if (Nf) *Nf = 4; if (v0) *v0 = quad_v0; if (J) *J = quad_J; if (invJ) *invJ = quad_invJ; if (detJ) *detJ = quad_detJ; break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unsupported polytope type %s", DMPolytopeTypes[ct]);
  }
  PetscFunctionReturn(0);
}

/* Gets the affine map from the original cell to each subcell */
static PetscErrorCode DMPlexCellRefinerGetAffineTransforms_Regular(DMPlexCellRefiner cr, DMPolytopeType ct, PetscInt *Nc, PetscReal *v0[], PetscReal *J[], PetscReal *invJ[])
{
  /*
   2
   |\
   | \
   |  \
   |   \
   | C  \
   |     \
   |      \
   2---1---1
   |\  D  / \
   | 2   0   \
   |A \ /  B  \
   0---0-------1
   */
  static PetscReal tri_v0[]   = {-1.0, -1.0,  0.0, -1.0,  -1.0, 0.0,  0.0, -1.0};
  static PetscReal tri_J[]    = {0.5, 0.0,
                                 0.0, 0.5,

                                 0.5, 0.0,
                                 0.0, 0.5,

                                 0.5, 0.0,
                                 0.0, 0.5,

                                 0.0, -0.5,
                                 0.5,  0.5};
  static PetscReal tri_invJ[] = {2.0, 0.0,
                                 0.0, 2.0,

                                 2.0, 0.0,
                                 0.0, 2.0,

                                 2.0, 0.0,
                                 0.0, 2.0,

                                 2.0,  2.0,
                                -2.0,  0.0};
    /*
     3---------2---------2
     |         |         |
     |    D    2    C    |
     |         |         |
     3----3----0----1----1
     |         |         |
     |    A    0    B    |
     |         |         |
     0---------0---------1
     */
  static PetscReal quad_v0[]   = {-1.0, -1.0,  0.0, -1.0,  0.0, 0.0,  -1.0, 0.0};
  static PetscReal quad_J[]    = {0.5, 0.0,
                                  0.0, 0.5,

                                  0.5, 0.0,
                                  0.0, 0.5,

                                  0.5, 0.0,
                                  0.0, 0.5,

                                  0.5, 0.0,
                                  0.0, 0.5};
  static PetscReal quad_invJ[] = {2.0, 0.0,
                                  0.0, 2.0,

                                  2.0, 0.0,
                                  0.0, 2.0,

                                  2.0, 0.0,
                                  0.0, 2.0,

                                  2.0, 0.0,
                                  0.0, 2.0};
    /*
     Bottom (viewed from top)    Top
     1---------2---------2       7---------2---------6
     |         |         |       |         |         |
     |    B    2    C    |       |    H    2    G    |
     |         |         |       |         |         |
     3----3----0----1----1       3----3----0----1----1
     |         |         |       |         |         |
     |    A    0    D    |       |    E    0    F    |
     |         |         |       |         |         |
     0---------0---------3       4---------0---------5
     */
  static PetscReal hex_v0[]   = {-1.0, -1.0, -1.0,  -1.0,  0.0, -1.0,  0.0, 0.0, -1.0,   0.0, -1.0, -1.0,
                                 -1.0, -1.0,  0.0,   0.0, -1.0,  0.0,  0.0, 0.0,  0.0,  -1.0,  0.0,  0.0};
  static PetscReal hex_J[]    = {0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5,

                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5,

                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5,

                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5,

                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5,

                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5,

                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5,

                                 0.5, 0.0, 0.0,
                                 0.0, 0.5, 0.0,
                                 0.0, 0.0, 0.5};
  static PetscReal hex_invJ[] = {2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0,

                                 2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0,

                                 2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0,

                                 2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0,

                                 2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0,

                                 2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0,

                                 2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0,

                                 2.0, 0.0, 0.0,
                                 0.0, 2.0, 0.0,
                                 0.0, 0.0, 2.0};

  PetscFunctionBegin;
  switch (ct) {
  case DM_POLYTOPE_TRIANGLE:      if (Nc) *Nc = 4; if (v0) *v0 = tri_v0;  if (J) *J = tri_J;  if (invJ) *invJ = tri_invJ;  break;
  case DM_POLYTOPE_QUADRILATERAL: if (Nc) *Nc = 4; if (v0) *v0 = quad_v0; if (J) *J = quad_J; if (invJ) *invJ = quad_invJ; break;
  case DM_POLYTOPE_HEXAHEDRON:    if (Nc) *Nc = 8; if (v0) *v0 = hex_v0;  if (J) *J = hex_J;  if (invJ) *invJ = hex_invJ;  break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unsupported polytope type %s", DMPolytopeTypes[ct]);
  }
  PetscFunctionReturn(0);
}

/* Should this be here or in the DualSpace somehow? */
PetscErrorCode CellRefinerInCellTest_Internal(DMPolytopeType ct, const PetscReal point[], PetscBool *inside)
{
  PetscReal sum = 0.0;
  PetscInt  d;

  PetscFunctionBegin;
  *inside = PETSC_TRUE;
  switch (ct) {
  case DM_POLYTOPE_TRIANGLE:
  case DM_POLYTOPE_TETRAHEDRON:
    for (d = 0; d < DMPolytopeTypeGetDim(ct); ++d) {
      if (point[d] < -1.0) {*inside = PETSC_FALSE; break;}
      sum += point[d];
    }
    if (sum > PETSC_SMALL) {*inside = PETSC_FALSE; break;}
    break;
  case DM_POLYTOPE_QUADRILATERAL:
  case DM_POLYTOPE_HEXAHEDRON:
    for (d = 0; d < DMPolytopeTypeGetDim(ct); ++d)
      if (PetscAbsReal(point[d]) > 1.+PETSC_SMALL) {*inside = PETSC_FALSE; break;}
    break;
  default:
    SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Unsupported polytope type %s", DMPolytopeTypes[ct]);
  }
  PetscFunctionReturn(0);
}

/* Regular Refinment of Hybrid Meshes

   We would like to express regular refinement as a small set of rules that can be applied on every point of the Plex
   to automatically generate a refined Plex. In fact, we would like these rules to be general enough to encompass other
   transformations, such as changing from one type of cell to another, as simplex to hex.

   To start, we can create a function that takes an original cell type and returns the number of new cells replacing it
   and the types of the new cells.

   We need the group multiplication table for group actions from the dihedral group for each cell type.

   We need an operator which takes in a cell, and produces a new set of cells with new faces and correct orientations. I think
   we can just write this operator for faces with identity, and then compose the face orientation actions to get the actual
   (face, orient) pairs for each subcell.
*/

/*
  Input Parameters:
+ ct - The type of the input cell
. co - The orientation of this cell in it parent
- cp - The requested cone point of this cell assuming orientation 0

  Output Parameters:
. cpnew - The new cone point, taking inout account the orientation co
*/
PETSC_STATIC_INLINE PetscErrorCode DMPolytopeMapCell(DMPolytopeType ct, PetscInt co, PetscInt cp, PetscInt *cpnew)
{
  const PetscInt csize = DMPolytopeTypeGetConeSize(ct);

  PetscFunctionBeginHot;
  if (ct == DM_POLYTOPE_POINT) {*cpnew = cp;}
  else                         {*cpnew = (co < 0 ? -(co+1)-cp+csize : co+cp)%csize;}
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerGetCellVertices_Regular(DMPlexCellRefiner cr, DMPolytopeType ct, PetscInt *Nv, PetscReal *subcellV[])
{
  static PetscReal seg_v[]  = {-1.0,  0.0,  1.0};
  static PetscReal tri_v[]  = {-1.0, -1.0,  1.0, -1.0,  -1.0, 1.0,  0.0, -1.0,  0.0, 0.0,  -1.0, 0.0};
  static PetscReal quad_v[] = {-1.0, -1.0,  1.0, -1.0,   1.0, 1.0,  -1.0, 1.0,  0.0, -1.0,  1.0, 0.0,   0.0, 1.0,  -1.0, 0.0,  0.0, 0.0};
  static PetscReal tet_v[]  = {-1.0, -1.0, -1.0,   0.0, -1.0, -1.0,   1.0, -1.0, -1.0,
                               -1.0,  0.0, -1.0,   0.0,  0.0, -1.0,  -1.0,  1.0, -1.0,
                               -1.0, -1.0,  0.0,   0.0, -1.0,  0.0,  -1.0,  0.0,  0.0,  -1.0, -1.0,  1.0};
  static PetscReal hex_v[]  = {-1.0, -1.0, -1.0,   0.0, -1.0, -1.0,   1.0, -1.0, -1.0,
                               -1.0,  0.0, -1.0,   0.0,  0.0, -1.0,   1.0,  0.0, -1.0,
                               -1.0,  1.0, -1.0,   0.0,  1.0, -1.0,   1.0,  1.0, -1.0,
                               -1.0, -1.0,  0.0,   0.0, -1.0,  0.0,   1.0, -1.0,  0.0,
                               -1.0,  0.0,  0.0,   0.0,  0.0,  0.0,   1.0,  0.0,  0.0,
                               -1.0,  1.0,  0.0,   0.0,  1.0,  0.0,   1.0,  1.0,  0.0,
                               -1.0, -1.0,  1.0,   0.0, -1.0,  1.0,   1.0, -1.0,  1.0,
                               -1.0,  0.0,  1.0,   0.0,  0.0,  1.0,   1.0,  0.0,  1.0,
                               -1.0,  1.0,  1.0,   0.0,  1.0,  1.0,   1.0,  1.0,  1.0};

  PetscFunctionBegin;
  switch (ct) {
    case DM_POLYTOPE_SEGMENT:       *Nv =  3; *subcellV = seg_v;  break;
    case DM_POLYTOPE_TRIANGLE:      *Nv =  6; *subcellV = tri_v;  break;
    case DM_POLYTOPE_QUADRILATERAL: *Nv =  9; *subcellV = quad_v; break;
    case DM_POLYTOPE_TETRAHEDRON:   *Nv = 10; *subcellV = tet_v;  break;
    case DM_POLYTOPE_HEXAHEDRON:    *Nv = 27; *subcellV = hex_v;  break;
    default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "No subcell vertices for cell type %s", DMPolytopeTypes[ct]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerGetCellVertices_ToBox(DMPlexCellRefiner cr, DMPolytopeType ct, PetscInt *Nv, PetscReal *subcellV[])
{
  static PetscReal tri_v[] = {-1.0, -1.0,  1.0, -1.0,  -1.0, 1.0,  0.0, -1.0,  0.0, 0.0,  -1.0, 0.0,  -1.0/3.0, -1.0/3.0};
  static PetscReal tet_v[] = {-1.0, -1.0, -1.0,   0.0, -1.0, -1.0,   1.0, -1.0, -1.0,
                              -1.0,  0.0, -1.0,  -1.0/3.0, -1.0/3.0, -1.0,   0.0,  0.0, -1.0,  -1.0,  1.0, -1.0,
                              -1.0, -1.0,  0.0,  -1.0/3.0, -1.0, -1.0/3.0,   0.0, -1.0,  0.0,
                              -1.0, -1.0/3.0, -1.0/3.0,  -1.0/3.0, -1.0/3.0, -1.0/3.0,  -1.0,  0.0,  0.0,
                              -1.0, -1.0,  1.0,  -0.5, -0.5, -0.5};
  PetscErrorCode   ierr;

  PetscFunctionBegin;
  switch (ct) {
    case DM_POLYTOPE_SEGMENT:
    case DM_POLYTOPE_QUADRILATERAL:
    case DM_POLYTOPE_HEXAHEDRON:
      ierr = DMPlexCellRefinerGetCellVertices_Regular(cr, ct, Nv, subcellV);CHKERRQ(ierr);break;
    case DM_POLYTOPE_TRIANGLE:    *Nv =  7; *subcellV = tri_v; break;
    case DM_POLYTOPE_TETRAHEDRON: *Nv = 15; *subcellV = tet_v;  break;
    default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "No subcell vertices for cell type %s", DMPolytopeTypes[ct]);
  }
  PetscFunctionReturn(0);
}

/*
  DMPlexCellRefinerGetCellVertices - Get the set of refined vertices lying in the closure of a reference cell of given type

  Input Parameters:
+ cr - The DMPlexCellRefiner object
- ct - The cell type

  Output Parameters:
+ Nv       - The number of refined vertices in the closure of the reference cell of given type
- subcellV - The coordinates of these vertices in the reference cell

  Level: developer

.seealso: DMPlexCellRefinerGetSubcellVertices()
*/
static PetscErrorCode DMPlexCellRefinerGetCellVertices(DMPlexCellRefiner cr, DMPolytopeType ct, PetscInt *Nv, PetscReal *subcellV[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*cr->ops->getcellvertices)(cr, ct, Nv, subcellV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerGetSubcellVertices_Regular(DMPlexCellRefiner cr, DMPolytopeType ct, DMPolytopeType rct, PetscInt r, PetscInt *Nv, PetscInt *subcellV[])
{
  static PetscInt seg_v[]  = {0, 1, 1, 2};
  static PetscInt tri_v[]  = {0, 3, 5,  3, 1, 4,  5, 4, 2,  3, 4, 5};
  static PetscInt quad_v[] = {0, 4, 8, 7,  4, 1, 5, 8,  8, 5, 2, 6,  7, 8, 6, 3};
  static PetscInt tet_v[]  = {0, 3, 1, 6,  3, 2, 4, 8,  1, 4, 5, 7,  6, 8, 7, 9,
                              1, 6, 3, 7,  8, 4, 3, 7,  7, 3, 1, 4,  7, 3, 8, 6};
  static PetscInt hex_v[]  = {0,  3,  4,  1,  9, 10, 13, 12,   3,  6,  7,  4, 12, 13, 16, 15,   4,  7,  8,  5, 13, 14, 17, 16,   1,  4 , 5 , 2, 10, 11, 14, 13,
                              9, 12, 13, 10, 18, 19, 22, 21,  10, 13, 14, 11, 19, 20, 23, 22,  13, 16, 17, 14, 22, 23, 26, 25,  12, 15, 16, 13, 21, 22, 25, 24};

  PetscFunctionBegin;
  if (ct != rct) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Cell type %s does not produce %s", DMPolytopeTypes[ct], DMPolytopeTypes[rct]);
  switch (ct) {
    case DM_POLYTOPE_SEGMENT:       *Nv = 2; *subcellV = &seg_v[r*(*Nv)];  break;
    case DM_POLYTOPE_TRIANGLE:      *Nv = 3; *subcellV = &tri_v[r*(*Nv)];  break;
    case DM_POLYTOPE_QUADRILATERAL: *Nv = 4; *subcellV = &quad_v[r*(*Nv)]; break;
    case DM_POLYTOPE_TETRAHEDRON:   *Nv = 4; *subcellV = &tet_v[r*(*Nv)];  break;
    case DM_POLYTOPE_HEXAHEDRON:    *Nv = 8; *subcellV = &hex_v[r*(*Nv)];  break;
    default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "No subcell vertices for cell type %s", DMPolytopeTypes[ct]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerGetSubcellVertices_ToBox(DMPlexCellRefiner cr, DMPolytopeType ct, DMPolytopeType rct, PetscInt r, PetscInt *Nv, PetscInt *subcellV[])
{
  static PetscInt tri_v[]  = {0, 3, 6, 5,  3, 1, 4, 6,  5, 6, 4, 2};
  static PetscInt tet_v[]  = {0,  3,  4,  1,  7,  8, 14, 10,   6, 12, 11,  5,  3,  4, 14, 10,   2,  5, 11,  9,  1,  8, 14,  4,  13, 12 , 10,  7,  9,  8, 14, 11};
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  switch (ct) {
    case DM_POLYTOPE_SEGMENT:
    case DM_POLYTOPE_QUADRILATERAL:
    case DM_POLYTOPE_HEXAHEDRON:
      ierr = DMPlexCellRefinerGetSubcellVertices_Regular(cr, ct, rct, r, Nv, subcellV);CHKERRQ(ierr);break;
    case DM_POLYTOPE_TRIANGLE:
      if (rct != DM_POLYTOPE_QUADRILATERAL) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Cell type %s does not produce %s", DMPolytopeTypes[ct], DMPolytopeTypes[rct]);
      *Nv = 4; *subcellV = &tri_v[r*(*Nv)]; break;
    case DM_POLYTOPE_TETRAHEDRON:
      if (rct != DM_POLYTOPE_HEXAHEDRON) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_INCOMP, "Cell type %s does not produce %s", DMPolytopeTypes[ct], DMPolytopeTypes[rct]);
      *Nv = 8; *subcellV = &tet_v[r*(*Nv)]; break;
    default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "No subcell vertices for cell type %s", DMPolytopeTypes[ct]);
  }
  PetscFunctionReturn(0);
}

/*
  DMPlexCellRefinerGetSubcellVertices - Get the set of refined vertices defining a subcell in the reference cell of given type

  Input Parameters:
+ cr  - The DMPlexCellRefiner object
. ct  - The cell type
. rct - The type of subcell
- r   - The subcell index

  Output Parameters:
+ Nv       - The number of refined vertices in the subcell
- subcellV - The indices of these vertices in the set of vertices returned by DMPlexCellRefinerGetCellVertices()

  Level: developer

.seealso: DMPlexCellRefinerGetCellVertices()
*/
static PetscErrorCode DMPlexCellRefinerGetSubcellVertices(DMPlexCellRefiner cr, DMPolytopeType ct, DMPolytopeType rct, PetscInt r, PetscInt *Nv, PetscInt *subcellV[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = (*cr->ops->getsubcellvertices)(cr, ct, rct, r, Nv, subcellV);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  Input Parameters:
+ cr  - The DMPlexCellRefiner
. pct - The cell type of the parent, from whom the new cell is being produced
. po  - The orientation of the parent cell in its enclosing parent
. ct  - The type being produced
. r   - The replica number requested for the produced cell type
- o   - The relative orientation of the replica

  Output Parameters:
+ rnew - The replica number, given the orientation of of the parent
- onew - The replica orientation, given the orientation of the parent
*/
static PetscErrorCode DMPlexCellRefinerMapSubcells(DMPlexCellRefiner cr, DMPolytopeType pct, PetscInt po, DMPolytopeType ct, PetscInt r, PetscInt o, PetscInt *rnew, PetscInt *onew)
{
  PetscErrorCode ierr;

  PetscFunctionBeginHot;
  ierr = (*cr->ops->mapsubcells)(cr, pct, po, ct, r, o, rnew, onew);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerMapSubcells_Regular(DMPlexCellRefiner cr, DMPolytopeType pct, PetscInt po, DMPolytopeType ct, PetscInt r, PetscInt o, PetscInt *rnew, PetscInt *onew)
{
  /* We shift any input orientation in order to make it non-negative
       The orientation array o[po][o] gives the orientation the new replica rnew has to have in order to reproduce the face sequence from (r, o)
       The replica array r[po][r] gives the new replica number rnew given that the parent point has orientation po
       Overall, replica (r, o) in a parent with orientation 0 matches replica (rnew, onew) in a parent with orientation po
  */
  PetscInt tri_seg_o[] = {-2, 0,
                          -2, 0,
                          -2, 0,
                          0, -2,
                          0, -2,
                          0, -2};
  PetscInt tri_seg_r[] = {1, 0, 2,
                          0, 2, 1,
                          2, 1, 0,
                          0, 1, 2,
                          1, 2, 0,
                          2, 0, 1};
  PetscInt tri_tri_o[] = {0,  1,  2, -3, -2, -1,
                          2,  0,  1, -2, -1, -3,
                          1,  2,  0, -1, -3, -2,
                         -3, -2, -1,  0,  1,  2,
                         -1, -3, -2,  1,  2,  0,
                         -2, -1, -3,  2,  0,  1};
  /* orientation if the replica is the center triangle */
  PetscInt tri_tri_o_c[] = {2,  0,  1, -2, -1, -3,
                            1,  2,  0, -1, -3, -2,
                            0,  1,  2, -3, -2, -1,
                           -3, -2, -1,  0,  1,  2,
                           -1, -3, -2,  1,  2,  0,
                           -2, -1, -3,  2,  0,  1};
  PetscInt tri_tri_r[] = {0, 2, 1, 3,
                          2, 1, 0, 3,
                          1, 0, 2, 3,
                          0, 1, 2, 3,
                          1, 2, 0, 3,
                          2, 0, 1, 3};
  PetscInt quad_seg_r[] = {3, 2, 1, 0,
                           2, 1, 0, 3,
                           1, 0, 3, 2,
                           0, 3, 2, 1,
                           0, 1, 2, 3,
                           1, 2, 3, 0,
                           2, 3, 0, 1,
                           3, 0, 1, 2};
  PetscInt quad_quad_o[] = { 0,  1,  2,  3, -4, -3, -2, -1,
                             4,  0,  1,  2, -3, -2, -1, -4,
                             3,  4,  0,  1, -2, -1, -4, -3,
                             2,  3,  4,  0, -1, -4, -3, -2,
                            -4, -3, -2, -1,  0,  1,  2,  3,
                            -1, -4, -3, -2,  1,  2,  3,  0,
                            -2, -1, -4, -3,  2,  3,  0,  1,
                            -3, -2, -1, -4,  3,  0,  1,  2};
  PetscInt quad_quad_r[] = {0, 3, 2, 1,
                            3, 2, 1, 0,
                            2, 1, 0, 3,
                            1, 0, 3, 2,
                            0, 1, 2, 3,
                            1, 2, 3, 0,
                            2, 3, 0, 1,
                            3, 0, 1, 2};
  PetscInt tquad_tquad_o[] = { 0,  1, -2, -1,
                               1,  0, -1, -2,
                              -2, -1,  0,  1,
                              -1, -2,  1,  0};
  PetscInt tquad_tquad_r[] = {1, 0,
                              1, 0,
                              0, 1,
                              0, 1};

  PetscFunctionBeginHot;
  /* The default is no transformation */
  *rnew = r;
  *onew = o;
  switch (pct) {
    case DM_POLYTOPE_SEGMENT:
      if (ct == DM_POLYTOPE_SEGMENT) {
        if      (po == 0 || po == -1) {*rnew = r;       *onew = o;}
        else if (po == 1 || po == -2) {*rnew = (r+1)%2; *onew = (o == 0 || o == -1) ? -2 : 0;}
        else SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Invalid orientation %D for segment", po);
      }
      break;
    case DM_POLYTOPE_TRIANGLE:
      switch (ct) {
        case DM_POLYTOPE_SEGMENT:
          if (o == -1) o = 0;
          if (o == -2) o = 1;
          *onew = tri_seg_o[(po+3)*2+o];
          *rnew = tri_seg_r[(po+3)*3+r];
          break;
        case DM_POLYTOPE_TRIANGLE:
          *onew = r == 3 && po < 0 ? tri_tri_o_c[((po+3)%3)*6+o+3] : tri_tri_o[(po+3)*6+o+3];
          *rnew = tri_tri_r[(po+3)*4+r];
          break;
        default: break;
      }
      break;
    case DM_POLYTOPE_QUADRILATERAL:
      switch (ct) {
        case DM_POLYTOPE_SEGMENT:
          *onew = o;
          *rnew = quad_seg_r[(po+4)*4+r];
          break;
        case DM_POLYTOPE_QUADRILATERAL:
          *onew = quad_quad_o[(po+4)*8+o+4];
          *rnew = quad_quad_r[(po+4)*4+r];
          break;
        default: break;
      }
      break;
    case DM_POLYTOPE_SEG_PRISM_TENSOR:
      switch (ct) {
        /* DM_POLYTOPE_POINT_PRISM_TENSOR does not change */
        case DM_POLYTOPE_SEG_PRISM_TENSOR:
          *onew = tquad_tquad_o[(po+2)*4+o+2];
          *rnew = tquad_tquad_r[(po+2)*2+r];
          break;
        default: break;
      }
      break;
    default: break;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerMapSubcells_ToBox(DMPlexCellRefiner cr, DMPolytopeType pct, PetscInt po, DMPolytopeType ct, PetscInt r, PetscInt o, PetscInt *rnew, PetscInt *onew)
{
  PetscErrorCode ierr;
  /* We shift any input orientation in order to make it non-negative
       The orientation array o[po][o] gives the orientation the new replica rnew has to have in order to reproduce the face sequence from (r, o)
       The replica array r[po][r] gives the new replica number rnew given that the parent point has orientation po
       Overall, replica (r, o) in a parent with orientation 0 matches replica (rnew, onew) in a parent with orientation po
  */
  PetscInt tri_seg_o[] = {0, -2,
                          0, -2,
                          0, -2,
                          0, -2,
                          0, -2,
                          0, -2};
  PetscInt tri_seg_r[] = {2, 1, 0,
                          1, 0, 2,
                          0, 2, 1,
                          0, 1, 2,
                          1, 2, 0,
                          2, 0, 1};
  PetscInt tri_tri_o[] = {0,  1,  2,  3, -4, -3, -2, -1,
                          3,  0,  1,  2, -3, -2, -1, -4,
                          1,  2,  3,  0, -1, -4, -3, -2,
                         -4, -3, -2, -1,  0,  1,  2,  3,
                         -1, -4, -3, -2,  1,  2,  3,  0,
                         -3, -2, -1, -4,  3,  0,  1,  2};
  PetscInt tri_tri_r[] = {0, 2, 1,
                          2, 1, 0,
                          1, 0, 2,
                          0, 1, 2,
                          1, 2, 0,
                          2, 0, 1};
  PetscInt tquad_tquad_o[] = { 0,  1, -2, -1,
                               1,  0, -1, -2,
                              -2, -1,  0,  1,
                              -1, -2,  1,  0};
  PetscInt tquad_tquad_r[] = {1, 0,
                              1, 0,
                              0, 1,
                              0, 1};
  PetscInt tquad_quad_o[] = {-2, -1, -4, -3,  2,  3,  0,  1,
                              1,  2,  3,  0, -1, -4, -3, -2,
                             -4, -3, -2, -1,  0,  1,  2,  3,
                              1,  0,  3,  2, -3, -4, -1, -2};

  PetscFunctionBeginHot;
  *rnew = r;
  *onew = o;
  switch (pct) {
    case DM_POLYTOPE_TRIANGLE:
      switch (ct) {
        case DM_POLYTOPE_SEGMENT:
          if (o == -1) o = 0;
          if (o == -2) o = 1;
          *onew = tri_seg_o[(po+3)*2+o];
          *rnew = tri_seg_r[(po+3)*3+r];
          break;
        case DM_POLYTOPE_QUADRILATERAL:
          *onew = tri_tri_o[(po+3)*8+o+4];
          /* TODO See sheet, for fixing po == 1 and 2 */
          if (po ==  2 && r == 2 && o >= 0) *onew = tri_tri_o[(po+3)*8+((o+3)%4)+4];
          if (po ==  2 && r == 2 && o <  0) *onew = tri_tri_o[(po+3)*8+((o+5)%4)];
          if (po ==  1 && r == 1 && o >= 0) *onew = tri_tri_o[(po+3)*8+((o+1)%4)+4];
          if (po ==  1 && r == 1 && o <  0) *onew = tri_tri_o[(po+3)*8+((o+7)%4)];
          if (po == -1 && r == 2 && o >= 0) *onew = tri_tri_o[(po+3)*8+((o+3)%4)+4];
          if (po == -1 && r == 2 && o <  0) *onew = tri_tri_o[(po+3)*8+((o+5)%4)];
          if (po == -2 && r == 1 && o >= 0) *onew = tri_tri_o[(po+3)*8+((o+1)%4)+4];
          if (po == -2 && r == 1 && o <  0) *onew = tri_tri_o[(po+3)*8+((o+7)%4)];
          *rnew = tri_tri_r[(po+3)*3+r];
          break;
        default: break;
      }
      break;
    case DM_POLYTOPE_SEG_PRISM_TENSOR:
      switch (ct) {
        /* DM_POLYTOPE_POINT_PRISM_TENSOR does not change */
        case DM_POLYTOPE_SEG_PRISM_TENSOR:
          *onew = tquad_tquad_o[(po+2)*4+o+2];
          *rnew = tquad_tquad_r[(po+2)*2+r];
          break;
        case DM_POLYTOPE_QUADRILATERAL:
          *onew = tquad_quad_o[(po+2)*8+o+4];
          *rnew = tquad_tquad_r[(po+2)*2+r];
          break;
        default: break;
      }
      break;
    default:
      ierr = DMPlexCellRefinerMapSubcells_Regular(cr, pct, po, ct, r, o, rnew, onew);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerMapSubcells_ToSimplex(DMPlexCellRefiner cr, DMPolytopeType pct, PetscInt po, DMPolytopeType ct, PetscInt r, PetscInt o, PetscInt *rnew, PetscInt *onew)
{
  return DMPlexCellRefinerMapSubcells_Regular(cr, pct, po, ct, r, o, rnew, onew);
}

/*@
  DMPlexCellRefinerRefine - Return a description of the refinement for a given cell type

  Input Parameter:
. source - The cell type for a source point

  Output Parameter:
+ Nt     - The number of cell types generated by refinement
. target - The cell types generated
. size   - The number of subcells of each type, ordered by dimension
. cone   - A list of the faces for each subcell of the same type as source
- ornt   - A list of the face orientations for each subcell of the same type as source

  Note: The cone array gives the cone of each subcell listed by the first three outputs. For the each cone point, we
  need the cell type, point identifier, and orientation within the subcell. The orientation is with respect to the canonical
  division (described in these outputs) of the cell in the original mesh. The point identifier is given by
$   the number of cones to be taken, or 0 for the current cell
$   the cell cone point number at each level from which it is subdivided
$   the replica number r of the subdivision.
  The orientation is with respect to the canonical cone orientation. For example, the prescription for edge division is
$   Nt     = 2
$   target = {DM_POLYTOPE_POINT, DM_POLYTOPE_SEGMENT}
$   size   = {1, 2}
$   cone   = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 0, 0,  DM_POLYTOPE_POINT, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0}
$   ornt   = {                         0,                       0,                        0,                          0}

  Level: developer

.seealso: DMPlexCellRefinerCreate(), DMPlexRefineUniform()
@*/
PetscErrorCode DMPlexCellRefinerRefine(DMPlexCellRefiner cr, DMPolytopeType source, PetscInt *Nt, DMPolytopeType *target[], PetscInt *size[], PetscInt *cone[], PetscInt *ornt[])
{
  PetscErrorCode ierr;

  PetscFunctionBeginHot;
  ierr = (*cr->ops->refine)(cr, source, Nt, target, size, cone, ornt);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerRefine_Regular(DMPlexCellRefiner cr, DMPolytopeType source, PetscInt *Nt, DMPolytopeType *target[], PetscInt *size[], PetscInt *cone[], PetscInt *ornt[])
{
  /* All vertices remain in the refined mesh */
  static DMPolytopeType vertexT[] = {DM_POLYTOPE_POINT};
  static PetscInt       vertexS[] = {1};
  static PetscInt       vertexC[] = {0};
  static PetscInt       vertexO[] = {0};
  /* Split all edges with a new vertex, making two new 2 edges
     0--0--0--1--1
  */
  static DMPolytopeType edgeT[]   = {DM_POLYTOPE_POINT, DM_POLYTOPE_SEGMENT};
  static PetscInt       edgeS[]   = {1, 2};
  static PetscInt       edgeC[]   = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 0, 0,  DM_POLYTOPE_POINT, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0};
  static PetscInt       edgeO[]   = {                         0,                       0,                        0,                          0};
  /* Do not split tensor edges */
  static DMPolytopeType tedgeT[]  = {DM_POLYTOPE_POINT_PRISM_TENSOR};
  static PetscInt       tedgeS[]  = {1};
  static PetscInt       tedgeC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0};
  static PetscInt       tedgeO[]  = {                         0,                          0};
  /* Add 3 edges inside every triangle, making 4 new triangles.
   2
   |\
   | \
   |  \
   0   1
   | C  \
   |     \
   |      \
   2---1---1
   |\  D  / \
   1 2   0   0
   |A \ /  B  \
   0-0-0---1---1
  */
  static DMPolytopeType triT[]    = {DM_POLYTOPE_SEGMENT, DM_POLYTOPE_TRIANGLE};
  static PetscInt       triS[]    = {3, 4};
  static PetscInt       triC[]    = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0,
                                     DM_POLYTOPE_POINT, 1, 1, 0, DM_POLYTOPE_POINT, 1, 2, 0,
                                     DM_POLYTOPE_POINT, 1, 2, 0, DM_POLYTOPE_POINT, 1, 0, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 2, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 0,    2};
  static PetscInt       triO[]    = {0, 0,
                                     0, 0,
                                     0, 0,
                                     0, -2,  0,
                                     0,  0, -2,
                                    -2,  0,  0,
                                     0,  0,  0};
  /* Add a vertex in the center of each quadrilateral, and 4 edges inside, making 4 new quads.
     3----1----2----0----2
     |         |         |
     0    D    2    C    1
     |         |         |
     3----3----0----1----1
     |         |         |
     1    A    0    B    0
     |         |         |
     0----0----0----1----1
  */
  static DMPolytopeType quadT[]   = {DM_POLYTOPE_POINT, DM_POLYTOPE_SEGMENT, DM_POLYTOPE_QUADRILATERAL};
  static PetscInt       quadS[]   = {1, 4, 4};
  static PetscInt       quadC[]   = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 1, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 2, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 3, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 1, 3, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 1, 2, 0, DM_POLYTOPE_SEGMENT, 0,    2,
                                     DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 2, 1, DM_POLYTOPE_SEGMENT, 1, 3, 0};
  static PetscInt       quadO[]   = {0, 0,
                                     0, 0,
                                     0, 0,
                                     0, 0,
                                     0,  0, -2,  0,
                                     0,  0,  0, -2,
                                    -2,  0,  0,  0,
                                     0, -2,  0,  0};
  /* Add 1 edge inside every tensor quad, making 2 new tensor quads
     2----2----1----3----3
     |         |         |
     |         |         |
     |         |         |
     4    A    6    B    5
     |         |         |
     |         |         |
     |         |         |
     0----0----0----1----1
  */
  static DMPolytopeType tquadT[]  = {DM_POLYTOPE_POINT_PRISM_TENSOR, DM_POLYTOPE_SEG_PRISM_TENSOR};
  static PetscInt       tquadS[]  = {1, 2};
  static PetscInt       tquadC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0,   0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_POINT_PRISM_TENSOR, 0,    0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 3, 0};
  static PetscInt       tquadO[]  = {0, 0,
                                     0, 0, 0, 0,
                                     0, 0, 0, 0};
  /* Add 1 edge and 8 triangles inside every cell, making 8 new tets
     The vertices of our reference tet are [(-1, -1, -1), (-1, 1, -1), (1, -1, -1), (-1, -1, 1)], which we call [v0, v1, v2, v3]. The first
     three edges are [v0, v1], [v1, v2], [v2, v0] called e0, e1, and e2, and then three edges to the top point [v0, v3], [v1, v3], [v2, v3]
     called e3, e4, and e5. The faces of a tet, given in DMPlexGetRawFaces_Internal() are
       [v0, v1, v2], [v0, v3, v1], [v0, v2, v3], [v2, v1, v3]
     The first four tets just cut off the corners, using the replica notation for new vertices,
       [v0,      (e0, 0), (e2, 0), (e3, 0)]
       [(e0, 0), v1,      (e1, 0), (e4, 0)]
       [(e2, 0), (e1, 0), v2,      (e5, 0)]
       [(e3, 0), (e4, 0), (e5, 0), v3     ]
     The next four tets match a vertex to the newly created faces from cutting off those first tets.
       [(e2, 0), (e3, 0), (e0, 0), (e5, 0)]
       [(e4, 0), (e1, 0), (e0, 0), (e5, 0)]
       [(e5, 0), (e0, 0), (e2, 0), (e1, 0)]
       [(e5, 0), (e0, 0), (e4, 0), (e3, 0)]
     We can see that a new edge is introduced in the cell [(e0, 0), (e5, 0)] which we call (-1, 0). The first four faces created are
       [(e2, 0), (e0, 0), (e3, 0)]
       [(e0, 0), (e1, 0), (e4, 0)]
       [(e2, 0), (e5, 0), (e1, 0)]
       [(e3, 0), (e4, 0), (e5, 0)]
     The next four, from the second group of tets, are
       [(e2, 0), (e0, 0), (e5, 0)]
       [(e4, 0), (e0, 0), (e5, 0)]
       [(e0, 0), (e1, 0), (e5, 0)]
       [(e5, 0), (e3, 0), (e0, 0)]
     I could write a program to generate these orientations by comparing the faces from GetRawFaces() with my existing table.
   */
  static DMPolytopeType tetT[]    = {DM_POLYTOPE_SEGMENT, DM_POLYTOPE_TRIANGLE, DM_POLYTOPE_TETRAHEDRON};
  static PetscInt       tetS[]    = {1, 8, 8};
  static PetscInt       tetC[]    = {DM_POLYTOPE_POINT, 2, 0, 0, 0, DM_POLYTOPE_POINT, 2, 2, 1, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_SEGMENT, 1, 2, 2,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 3, 0, DM_POLYTOPE_SEGMENT, 1, 1, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 0, DM_POLYTOPE_SEGMENT, 1, 3, 2, DM_POLYTOPE_SEGMENT, 1, 0, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 1, 3, 1, DM_POLYTOPE_SEGMENT, 1, 2, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 3, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 3, 2, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 1, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 0, DM_POLYTOPE_TRIANGLE, 1, 1, 0, DM_POLYTOPE_TRIANGLE, 1, 2, 0, DM_POLYTOPE_TRIANGLE, 0,    0,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 1, DM_POLYTOPE_TRIANGLE, 1, 1, 2, DM_POLYTOPE_TRIANGLE, 0,    1, DM_POLYTOPE_TRIANGLE, 1, 3, 1,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 2, DM_POLYTOPE_TRIANGLE, 0,    2, DM_POLYTOPE_TRIANGLE, 1, 2, 1, DM_POLYTOPE_TRIANGLE, 1, 3, 0,
                                     DM_POLYTOPE_TRIANGLE, 0,    3, DM_POLYTOPE_TRIANGLE, 1, 1, 1, DM_POLYTOPE_TRIANGLE, 1, 2, 2, DM_POLYTOPE_TRIANGLE, 1, 3, 2,
                                     DM_POLYTOPE_TRIANGLE, 0,    0, DM_POLYTOPE_TRIANGLE, 1, 2, 3, DM_POLYTOPE_TRIANGLE, 0,    4, DM_POLYTOPE_TRIANGLE, 0,    7,
                                     DM_POLYTOPE_TRIANGLE, 0,    1, DM_POLYTOPE_TRIANGLE, 1, 3, 3, DM_POLYTOPE_TRIANGLE, 0,    5, DM_POLYTOPE_TRIANGLE, 0,    6,
                                     DM_POLYTOPE_TRIANGLE, 0,    4, DM_POLYTOPE_TRIANGLE, 0,    6, DM_POLYTOPE_TRIANGLE, 0,    2, DM_POLYTOPE_TRIANGLE, 1, 0, 3,
                                     DM_POLYTOPE_TRIANGLE, 0,    5, DM_POLYTOPE_TRIANGLE, 0,    7, DM_POLYTOPE_TRIANGLE, 0,    3, DM_POLYTOPE_TRIANGLE, 1, 1, 3};
  static PetscInt       tetO[]    = {0, 0,
                                     0,  0,  0,
                                     0,  0,  0,
                                     0,  0,  0,
                                     0,  0,  0,
                                     0,  0, -2,
                                     0,  0, -2,
                                     0, -2, -2,
                                     0, -2,  0,
                                     0,  0,  0,  0,
                                     0,  0,  0,  0,
                                     0,  0,  0,  0,
                                     0,  0,  0,  0,
                                    -3,  0,  0, -2,
                                    -2,  1,  0,  0,
                                    -2, -2, -1,  2,
                                    -2,  0, -2,  1};
  /* Add a vertex in the center of each cell, add 6 edges and 12 quads inside every cell, making 8 new hexes
     The vertices of our reference hex are (-1, -1, -1), (-1, 1, -1), (1, 1, -1), (1, -1, -1), (-1, -1, 1), (1, -1, 1), (1, 1, 1), (-1, 1, 1) which we call [v0, v1, v2, v3, v4, v5, v6, v7]. The fours edges around the bottom [v0, v1], [v1, v2], [v2, v3], [v3, v0] are [e0, e1, e2, e3], and likewise around the top [v4, v5], [v5, v6], [v6, v7], [v7, v4] are [e4, e5, e6, e7]. Finally [v0, v4], [v1, v7], [v2, v6], [v3, v5] are [e9, e10, e11, e8]. The faces of a hex, given in DMPlexGetRawFaces_Internal(), oriented with outward normals, are
       [v0, v1, v2, v3] f0 bottom
       [v4, v5, v6, v7] f1 top
       [v0, v3, v5, v4] f2 front
       [v2, v1, v7, v6] f3 back
       [v3, v2, v6, v5] f4 right
       [v0, v4, v7, v1] f5 left
     The eight hexes are divided into four on the bottom, and four on the top,
       [v0,      (e0, 0),  (f0, 0),  (e3, 0),  (e9, 0), (f2, 0),  (c0, 0),  (f5, 0)]
       [(e0, 0), v1,       (e1, 0),  (f0, 0),  (f5, 0), (c0, 0),  (f3, 0),  (e10, 0)]
       [(f0, 0), (e1, 0),  v2,       (e2, 0),  (c0, 0), (f4, 0),  (e11, 0), (f3, 0)]
       [(e3, 0), (f0, 0),  (e2, 0),  v3,       (f2, 0), (e8, 0),  (f4, 0),  (c0, 0)]
       [(e9, 0), (f5, 0),  (c0, 0),  (f2, 0),  v4,      (e4, 0),  (f1, 0),  (e7, 0)]
       [(f2, 0), (c0, 0),  (f4, 0),  (e8, 0),  (e4, 0), v5,       (e5, 0),  (f1, 0)]
       [(c0, 0), (f3, 0),  (e11, 0), (f4, 0),  (f1, 0), (e5, 0),  v6,       (e6, 0)]
       [(f5, 0), (e10, 0), (f3, 0),  (c0, 0),  (e7, 0), (f1, 0),  (e6, 0),  v7]
     The 6 internal edges will go from the faces to the central vertex. The 12 internal faces can be divided into groups of 4 by the plane on which they sit. First the faces on the x-y plane are,
       [(e9, 0), (f2, 0),  (c0, 0),  (f5, 0)]
       [(f5, 0), (c0, 0),  (f3, 0),  (e10, 0)]
       [(c0, 0), (f4, 0),  (e11, 0), (f3, 0)]
       [(f2, 0), (e8, 0),  (f4, 0),  (c0, 0)]
     and on the x-z plane,
       [(f0, 0), (e0, 0), (f5, 0), (c0, 0)]
       [(c0, 0), (f5, 0), (e7, 0), (f1, 0)]
       [(f4, 0), (c0, 0), (f1, 0), (e5, 0)]
       [(e2, 0), (f0, 0), (c0, 0), (f4, 0)]
     and on the y-z plane,
       [(e3, 0), (f2, 0), (c0, 0), (f0, 0)]
       [(f2, 0), (e4, 0), (f1, 0), (c0, 0)]
       [(c0, 0), (f1, 0), (e6, 0), (f3, 0)]
       [(f0, 0), (c0, 0), (f3, 0), (e1, 0)]
  */
  static DMPolytopeType hexT[]    = {DM_POLYTOPE_POINT, DM_POLYTOPE_SEGMENT, DM_POLYTOPE_QUADRILATERAL, DM_POLYTOPE_HEXAHEDRON};
  static PetscInt       hexS[]    = {1, 6, 12, 8};
  static PetscInt       hexC[]    = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 1, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 2, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 3, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 4, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 5, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 3, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 0,    5, DM_POLYTOPE_SEGMENT, 1, 5, 0,
                                     DM_POLYTOPE_SEGMENT, 0,    5, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 1, 3, 1, DM_POLYTOPE_SEGMENT, 1, 5, 2,
                                     DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 1, 4, 1, DM_POLYTOPE_SEGMENT, 1, 3, 3, DM_POLYTOPE_SEGMENT, 0,    3,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 1, DM_POLYTOPE_SEGMENT, 1, 4, 3, DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 0,    2,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 5, 3, DM_POLYTOPE_SEGMENT, 0,    5, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 0,    5, DM_POLYTOPE_SEGMENT, 1, 5, 1, DM_POLYTOPE_SEGMENT, 1, 1, 3, DM_POLYTOPE_SEGMENT, 0,    1,
                                     DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 1, 4, 2,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 1, 4, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 0, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 0, 3,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 2, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 0,    2,
                                     DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_SEGMENT, 1, 3, 2, DM_POLYTOPE_SEGMENT, 0,    3,
                                     DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 1, 3, 0, DM_POLYTOPE_SEGMENT, 1, 0, 1,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 0, DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 0, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 0,    8, DM_POLYTOPE_QUADRILATERAL, 1, 5, 0,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 1, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 1, 3, 1, DM_POLYTOPE_QUADRILATERAL, 0,   11, DM_POLYTOPE_QUADRILATERAL, 1, 5, 3,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 2, DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 0,    7, DM_POLYTOPE_QUADRILATERAL, 1, 3, 0, DM_POLYTOPE_QUADRILATERAL, 1, 4, 1, DM_POLYTOPE_QUADRILATERAL, 0,   11,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 3, DM_POLYTOPE_QUADRILATERAL, 0,    3, DM_POLYTOPE_QUADRILATERAL, 1, 2, 1, DM_POLYTOPE_QUADRILATERAL, 0,    7, DM_POLYTOPE_QUADRILATERAL, 1, 4, 0, DM_POLYTOPE_QUADRILATERAL, 0,    8,
                                     DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 1, 0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 3, DM_POLYTOPE_QUADRILATERAL, 0,    5, DM_POLYTOPE_QUADRILATERAL, 0,    9, DM_POLYTOPE_QUADRILATERAL, 1, 5, 1,
                                     DM_POLYTOPE_QUADRILATERAL, 0,    3, DM_POLYTOPE_QUADRILATERAL, 1, 1, 1, DM_POLYTOPE_QUADRILATERAL, 1, 2, 2, DM_POLYTOPE_QUADRILATERAL, 0,    6, DM_POLYTOPE_QUADRILATERAL, 1, 4, 3, DM_POLYTOPE_QUADRILATERAL, 0,    9,
                                     DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 1, 1, 2, DM_POLYTOPE_QUADRILATERAL, 0,    6, DM_POLYTOPE_QUADRILATERAL, 1, 3, 3, DM_POLYTOPE_QUADRILATERAL, 1, 4, 2, DM_POLYTOPE_QUADRILATERAL, 0,   10,
                                     DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1, 1, 3, DM_POLYTOPE_QUADRILATERAL, 0,    5, DM_POLYTOPE_QUADRILATERAL, 1, 3, 2, DM_POLYTOPE_QUADRILATERAL, 0,   10, DM_POLYTOPE_QUADRILATERAL, 1, 5, 2};
  static PetscInt       hexO[]    = {0, 0,
                                     0, 0,
                                     0, 0,
                                     0, 0,
                                     0, 0,
                                     0, 0,
                                     0,  0, -2, -2,
                                     0, -2, -2,  0,
                                    -2, -2,  0,  0,
                                    -2,  0,  0, -2,
                                    -2,  0,  0, -2,
                                    -2, -2,  0,  0,
                                     0, -2, -2,  0,
                                     0,  0, -2, -2,
                                     0,  0, -2, -2,
                                    -2,  0,  0, -2,
                                    -2, -2,  0,  0,
                                     0, -2, -2,  0,
                                     0, 0,  0, 0, -4, 0,
                                     0, 0, -1, 0, -4, 0,
                                     0, 0, -1, 0,  0, 0,
                                     0, 0,  0, 0,  0, 0,
                                    -4, 0,  0, 0, -4, 0,
                                    -4, 0,  0, 0,  0, 0,
                                    -4, 0, -1, 0,  0, 0,
                                    -4, 0, -1, 0, -4, 0};
  /* Add 3 quads inside every triangular prism, making 4 new prisms. */
  static DMPolytopeType tripT[]   = {DM_POLYTOPE_SEGMENT, DM_POLYTOPE_TRIANGLE, DM_POLYTOPE_QUADRILATERAL, DM_POLYTOPE_TRI_PRISM};
  static PetscInt       tripS[]   = {3, 4, 6, 8};
  static PetscInt       tripC[]   = {DM_POLYTOPE_POINT, 1, 2, 0, DM_POLYTOPE_POINT, 1, 3, 0,
                                     DM_POLYTOPE_POINT, 1, 3, 0, DM_POLYTOPE_POINT, 1, 4, 0,
                                     DM_POLYTOPE_POINT, 1, 4, 0, DM_POLYTOPE_POINT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 3, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 4, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 1, DM_POLYTOPE_SEGMENT, 1, 3, 3, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 3, 1, DM_POLYTOPE_SEGMENT, 1, 4, 3,
                                     DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 0,    2,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 3, 0, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 4, 0, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 3, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 1, 2, 0, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 4, 0,
                                     DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 3, 2, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 1, 2, 2,
                                     DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 4, 2, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 1, 3, 2,
                                     DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 2, 2, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_SEGMENT, 1, 4, 2,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 0, DM_POLYTOPE_TRIANGLE, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 0, DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 1, 4, 1,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 2, DM_POLYTOPE_TRIANGLE, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1, 2, 1, DM_POLYTOPE_QUADRILATERAL, 1, 3, 0, DM_POLYTOPE_QUADRILATERAL, 0,    0,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 1, DM_POLYTOPE_TRIANGLE, 0,    2, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1, 3, 1, DM_POLYTOPE_QUADRILATERAL, 1, 4, 0,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 3, DM_POLYTOPE_TRIANGLE, 0,    3, DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 0,    2,
                                     DM_POLYTOPE_TRIANGLE, 0,    0, DM_POLYTOPE_TRIANGLE, 1, 1, 0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 3, DM_POLYTOPE_QUADRILATERAL, 0,    5, DM_POLYTOPE_QUADRILATERAL, 1, 4, 2,
                                     DM_POLYTOPE_TRIANGLE, 0,    1, DM_POLYTOPE_TRIANGLE, 1, 1, 1, DM_POLYTOPE_QUADRILATERAL, 1, 2, 2, DM_POLYTOPE_QUADRILATERAL, 1, 3, 3, DM_POLYTOPE_QUADRILATERAL, 0,    3,
                                     DM_POLYTOPE_TRIANGLE, 0,    2, DM_POLYTOPE_TRIANGLE, 1, 1, 2, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 1, 3, 2, DM_POLYTOPE_QUADRILATERAL, 1, 4, 3,
                                     DM_POLYTOPE_TRIANGLE, 0,    3, DM_POLYTOPE_TRIANGLE, 1, 1, 3, DM_POLYTOPE_QUADRILATERAL, 0,    3, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 0,    5};
  static PetscInt       tripO[]   = {0, 0,
                                     0, 0,
                                     0, 0,
                                     0, -2, -2,
                                    -2,  0, -2,
                                    -2, -2,  0,
                                     0,  0,  0,
                                    -2,  0, -2, -2,
                                    -2,  0, -2, -2,
                                    -2,  0, -2, -2,
                                     0, -2, -2,  0,
                                     0, -2, -2,  0,
                                     0, -2, -2,  0,
                                     0,  0,  0, -1,  0,
                                     0,  0,  0,  0, -1,
                                     0,  0, -1,  0,  0,
                                     2,  0,  0,  0,  0,
                                    -3,  0,  0, -1,  0,
                                    -3,  0,  0,  0, -1,
                                    -3,  0, -1,  0,  0,
                                    -3,  0,  0,  0,  0};
  /* Add 3 tensor quads inside every tensor triangular prism, making 4 new prisms.
      2
      |\
      | \
      |  \
      0---1

      2

      0   1

      2
      |\
      | \
      |  \
      0---1
  */
  static DMPolytopeType ttripT[]  = {DM_POLYTOPE_SEG_PRISM_TENSOR, DM_POLYTOPE_TRI_PRISM_TENSOR};
  static PetscInt       ttripS[]  = {3, 4};
  static PetscInt       ttripC[]  = {DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 3, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 4, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 4, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 2, 0,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 0, DM_POLYTOPE_TRIANGLE, 1, 1, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 1,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 1, DM_POLYTOPE_TRIANGLE, 1, 1, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 2, DM_POLYTOPE_TRIANGLE, 1, 1, 2, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 0,
                                     DM_POLYTOPE_TRIANGLE, 1, 0, 3, DM_POLYTOPE_TRIANGLE, 1, 1, 3, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,     1, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2};
  static PetscInt       ttripO[]  = {0, 0, 0, 0,
                                     0, 0, 0, 0,
                                     0, 0, 0, 0,
                                     0, 0,  0, -1,  0,
                                     0, 0,  0,  0, -1,
                                     0, 0, -1,  0,  0,
                                     0, 0,  0,  0,  0};
  /* Add 1 edge and 4 tensor quads inside every tensor quad prism, making 4 new prisms. */
  static DMPolytopeType tquadpT[]  = {DM_POLYTOPE_POINT_PRISM_TENSOR, DM_POLYTOPE_SEG_PRISM_TENSOR, DM_POLYTOPE_QUAD_PRISM_TENSOR};
  static PetscInt       tquadpS[]  = {1, 4, 4};
  static PetscInt       tquadpC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 4, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 3, DM_POLYTOPE_SEGMENT, 1, 1, 3, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 5, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 0, DM_POLYTOPE_QUADRILATERAL, 1, 1, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    3, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 5, 1,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 1, DM_POLYTOPE_QUADRILATERAL, 1, 1, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    1, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 2, DM_POLYTOPE_QUADRILATERAL, 1, 1, 2, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 3, DM_POLYTOPE_QUADRILATERAL, 1, 1, 3, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    3, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 5, 0};
  static PetscInt       tquadpO[]  = {0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0,  0,  0, -1,  0,
                                      0, 0,  0,  0,  0, -1,
                                      0, 0, -1,  0,  0,  0,
                                      0, 0,  0, -1,  0,  0};

  PetscFunctionBegin;
  switch (source) {
    case DM_POLYTOPE_POINT:              *Nt = 1; *target = vertexT; *size = vertexS; *cone = vertexC; *ornt = vertexO; break;
    case DM_POLYTOPE_SEGMENT:            *Nt = 2; *target = edgeT;   *size = edgeS;   *cone = edgeC;   *ornt = edgeO;   break;
    case DM_POLYTOPE_POINT_PRISM_TENSOR: *Nt = 1; *target = tedgeT;  *size = tedgeS;  *cone = tedgeC;  *ornt = tedgeO;  break;
    case DM_POLYTOPE_TRIANGLE:           *Nt = 2; *target = triT;    *size = triS;    *cone = triC;    *ornt = triO;    break;
    case DM_POLYTOPE_QUADRILATERAL:      *Nt = 3; *target = quadT;   *size = quadS;   *cone = quadC;   *ornt = quadO;   break;
    case DM_POLYTOPE_SEG_PRISM_TENSOR:   *Nt = 2; *target = tquadT;  *size = tquadS;  *cone = tquadC;  *ornt = tquadO;  break;
    case DM_POLYTOPE_TETRAHEDRON:        *Nt = 3; *target = tetT;    *size = tetS;    *cone = tetC;    *ornt = tetO;    break;
    case DM_POLYTOPE_HEXAHEDRON:         *Nt = 4; *target = hexT;    *size = hexS;    *cone = hexC;    *ornt = hexO;    break;
    case DM_POLYTOPE_TRI_PRISM:          *Nt = 4; *target = tripT;   *size = tripS;   *cone = tripC;   *ornt = tripO;   break;
    case DM_POLYTOPE_TRI_PRISM_TENSOR:   *Nt = 2; *target = ttripT;  *size = ttripS;  *cone = ttripC;  *ornt = ttripO;  break;
    case DM_POLYTOPE_QUAD_PRISM_TENSOR:  *Nt = 3; *target = tquadpT; *size = tquadpS; *cone = tquadpC; *ornt = tquadpO; break;
    default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "No refinement strategy for %s", DMPolytopeTypes[source]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerRefine_ToBox(DMPlexCellRefiner cr, DMPolytopeType source, PetscInt *Nt, DMPolytopeType *target[], PetscInt *size[], PetscInt *cone[], PetscInt *ornt[])
{
  PetscErrorCode ierr;
  /* Change tensor edges to segments */
  static DMPolytopeType tedgeT[]  = {DM_POLYTOPE_SEGMENT};
  static PetscInt       tedgeS[]  = {1};
  static PetscInt       tedgeC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0};
  static PetscInt       tedgeO[]  = {                         0,                          0};
  /* Add 1 vertex, 3 edges inside every triangle, making 3 new quadrilaterals.
   2
   |\
   | \
   |  \
   |   \
   0    1
   |     \
   |      \
   2       1
   |\     / \
   | 2   1   \
   |  \ /     \
   1   |       0
   |   0        \
   |   |         \
   |   |          \
   0-0-0-----1-----1
  */
  static DMPolytopeType triT[]    = {DM_POLYTOPE_POINT, DM_POLYTOPE_SEGMENT, DM_POLYTOPE_QUADRILATERAL};
  static PetscInt       triS[]    = {1, 3, 3};
  static PetscInt       triC[]    = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 0,    0,
                                     DM_POLYTOPE_POINT, 1, 1, 0, DM_POLYTOPE_POINT, 0,    0,
                                     DM_POLYTOPE_POINT, 1, 2, 0, DM_POLYTOPE_POINT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 2, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 1, 2, 0};
  static PetscInt       triO[]    = {0, 0,
                                     0, 0,
                                     0, 0,
                                     0,  0, -2,  0,
                                     0,  0,  0, -2,
                                     0, -2,  0,  0};
  /* Add 1 edge inside every tensor quad, making 2 new quadrilaterals
     2----2----1----3----3
     |         |         |
     |         |         |
     |         |         |
     4    A    6    B    5
     |         |         |
     |         |         |
     |         |         |
     0----0----0----1----1
  */
  static DMPolytopeType tquadT[]  = {DM_POLYTOPE_SEGMENT, DM_POLYTOPE_QUADRILATERAL};
  static PetscInt       tquadS[]  = {1, 2};
  static PetscInt       tquadC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0,
                                     /* TODO  Fix these */
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 3, 0, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 0,    0};
  static PetscInt       tquadO[]  = {0, 0,
                                     0, 0, -2, -2,
                                     0, 0, -2, -2};
  /* Add 6 triangles inside every cell, making 4 new hexs
     TODO: Need different SubcellMap(). Need to make a struct with the function pointers in it
     The vertices of our reference tet are [(-1, -1, -1), (-1, 1, -1), (1, -1, -1), (-1, -1, 1)], which we call [v0, v1, v2, v3]. The first
     three edges are [v0, v1], [v1, v2], [v2, v0] called e0, e1, and e2, and then three edges to the top point [v0, v3], [v1, v3], [v2, v3]
     called e3, e4, and e5. The faces of a tet, given in DMPlexGetRawFaces_Internal() are
       [v0, v1, v2], [v0, v3, v1], [v0, v2, v3], [v2, v1, v3]
     We make a new hex in each corner
       [v0, (e0, 0), (f0, 0), (e2, 0), (e3, 0), (f2, 0), (c0, 0), (f1, 0)]
       [v1, (e4, 0), (f3, 0), (e1, 0), (e0, 0), (f0, 0), (c0, 0), (f1, 0)]
       [v2, (e1, 0), (f3, 0), (e5, 0), (e2, 0), (f2, 0), (c0, 0), (f0, 0)]
       [v3, (e4, 0), (f1, 0), (e3, 0), (e5, 0), (f2, 0), (c0, 0), (f3, 0)]
     We create a new face for each edge
       [(e3, 0), (f2, 0), (c0, 0), (f1, 0)]
       [(f0, 0), (e0, 0), (f1, 0), (c0, 0)]
       [(e2, 0), (f0, 0), (c0, 0), (f2, 0)]
       [(f3, 0), (e4, 0), (f1, 0), (c0, 0)]
       [(e1, 0), (f3, 0), (c0, 0), (f0, 0)]
       [(e5, 0), (f3, 0), (c0, 0), (f2, 0)]
     I could write a program to generate these from the first hex by acting with the symmetry group to take one subcell into another.
   */
  static DMPolytopeType tetT[]    = {DM_POLYTOPE_POINT, DM_POLYTOPE_SEGMENT, DM_POLYTOPE_QUADRILATERAL, DM_POLYTOPE_HEXAHEDRON};
  static PetscInt       tetS[]    = {1, 4, 6, 4};
  static PetscInt       tetC[]    = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 1, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 2, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 3, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 2, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 0,    0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 3, 1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 0,    3,
                                     DM_POLYTOPE_SEGMENT, 1, 3, 0, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 0, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 3, 2, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 2, 1,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 0, DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 0, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 1, 1, 0,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 3, 1, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1, 0, 1, DM_POLYTOPE_QUADRILATERAL, 0,    3, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 1, 1, 2,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 3, 0, DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 1, 2, 1, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 0,    5, DM_POLYTOPE_QUADRILATERAL, 1, 0, 2,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 1, 1, DM_POLYTOPE_QUADRILATERAL, 0,    5, DM_POLYTOPE_QUADRILATERAL, 1, 2, 2, DM_POLYTOPE_QUADRILATERAL, 0,    3, DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 3, 2};
  static PetscInt       tetO[]    = {0, 0,
                                     0, 0,
                                     0, 0,
                                     0, 0,
                                     0,  0, -2, -2,
                                    -2,  0,  0, -2,
                                     0,  0, -2, -2,
                                    -2,  0,  0, -2,
                                     0,  0, -2, -2,
                                     0,  0, -2, -2,
                                     0,  0,  0,  0,  0,  0,
                                     1, -1,  1,  0,  0,  3,
                                     0, -4,  1, -1,  0,  3,
                                     1, -4,  3, -2, -4,  3};
  /* Add 3 quads inside every triangular prism, making 4 new prisms. */
  static DMPolytopeType tripT[]   = {DM_POLYTOPE_POINT, DM_POLYTOPE_SEGMENT, DM_POLYTOPE_QUADRILATERAL, DM_POLYTOPE_HEXAHEDRON};
  static PetscInt       tripS[]   = {1, 5, 9, 6};
  static PetscInt       tripC[]   = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 1, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 2, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 3, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_POINT, 1, 4, 0, DM_POLYTOPE_POINT, 0, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 3, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 1, 4, 1,
                                     DM_POLYTOPE_SEGMENT, 1, 2, 1, DM_POLYTOPE_SEGMENT, 1, 3, 3, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 0,    2,
                                     DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 1, 3, 1, DM_POLYTOPE_SEGMENT, 1, 4, 3,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 1, 3, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 1, 4, 0,
                                     DM_POLYTOPE_SEGMENT, 0,    2, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 1, 2, 2,
                                     DM_POLYTOPE_SEGMENT, 0,    3, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 1, 3, 2,
                                     DM_POLYTOPE_SEGMENT, 0,    4, DM_POLYTOPE_SEGMENT, 0,    1, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_SEGMENT, 1, 4, 2,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 0, DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 0, DM_POLYTOPE_QUADRILATERAL, 0,    5, DM_POLYTOPE_QUADRILATERAL, 0,    3, DM_POLYTOPE_QUADRILATERAL, 1, 4, 1,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 2, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1, 2, 1, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 1, 3, 0, DM_POLYTOPE_QUADRILATERAL, 0,    3,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 1, DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 0,    5, DM_POLYTOPE_QUADRILATERAL, 1, 3, 1, DM_POLYTOPE_QUADRILATERAL, 0,    4, DM_POLYTOPE_QUADRILATERAL, 1, 4, 0,
                                     DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 1, 0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 3, DM_POLYTOPE_QUADRILATERAL, 0,    8, DM_POLYTOPE_QUADRILATERAL, 0,    6, DM_POLYTOPE_QUADRILATERAL, 1, 4, 2,
                                     DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1, 1, 1, DM_POLYTOPE_QUADRILATERAL, 1, 2, 2, DM_POLYTOPE_QUADRILATERAL, 0,    7, DM_POLYTOPE_QUADRILATERAL, 1, 3, 3, DM_POLYTOPE_QUADRILATERAL, 0,    6,
                                     DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 1, 1, 2, DM_POLYTOPE_QUADRILATERAL, 0,    8, DM_POLYTOPE_QUADRILATERAL, 1, 3, 2, DM_POLYTOPE_QUADRILATERAL, 0,    7, DM_POLYTOPE_QUADRILATERAL, 1, 4, 3};
  static PetscInt       tripO[]   = {0, 0,
                                     0, 0,
                                     0, 0,
                                     0, 0,
                                     0, 0,
                                     0,  0, -2, -2,
                                    -2,  0,  0, -2,
                                     0, -2, -2,  0,
                                     0,  0, -2, -2,
                                     0,  0, -2, -2,
                                     0,  0, -2, -2,
                                     0, -2, -2,  0,
                                     0, -2, -2,  0,
                                     0, -2, -2,  0,
                                     0,  0,  0, -1,  0,  1,
                                     0,  0,  0,  0,  0, -4,
                                     0,  0,  0,  0, -1,  1,
                                    -4,  0,  0, -1,  0,  1,
                                    -4,  0,  0,  0,  0, -4,
                                    -4,  0,  0,  0, -1,  1};
  /* Add 3 tensor quads inside every tensor triangular prism, making 4 new tensor triangular prisms.
      2
      |\
      | \
      |  \
      0---1

      2

      0   1

      2
      |\
      | \
      |  \
      0---1
  */
  static DMPolytopeType ttripT[]  = {DM_POLYTOPE_POINT_PRISM_TENSOR, DM_POLYTOPE_SEG_PRISM_TENSOR, DM_POLYTOPE_QUAD_PRISM_TENSOR};
  static PetscInt       ttripS[]  = {1, 3, 3};
  static PetscInt       ttripC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 4, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 0, DM_POLYTOPE_QUADRILATERAL, 1, 1, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 1,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 1, DM_POLYTOPE_QUADRILATERAL, 1, 1, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    1, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 2, DM_POLYTOPE_QUADRILATERAL, 1, 1, 2, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 0};
  static PetscInt       ttripO[]  = {0, 0,
                                     0, 0, 0, 0,
                                     0, 0, 0, 0,
                                     0, 0, 0, 0,
                                     0, 0, 0,  0, -1, 0,
                                     0, 0, 0,  0,  0, -1,
                                     0, 0, 0, -1,  0, 0};
  /* TODO Add 3 quads inside every tensor triangular prism, making 4 new triangular prisms.
      2
      |\
      | \
      |  \
      0---1

      2

      0   1

      2
      |\
      | \
      |  \
      0---1
  */
  static DMPolytopeType ctripT[]  = {DM_POLYTOPE_SEGMENT, DM_POLYTOPE_QUADRILATERAL, DM_POLYTOPE_HEXAHEDRON};
  static PetscInt       ctripS[]  = {1, 3, 3};
  static PetscInt       ctripC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_SEGMENT, 1, 2, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_SEGMENT, 1, 3, 0,
                                     DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 0,    0, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_SEGMENT, 1, 4, 0,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 0, DM_POLYTOPE_QUADRILATERAL, 1, 1, 0, DM_POLYTOPE_QUADRILATERAL, 1, 2, 0, DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 0,    0, DM_POLYTOPE_QUADRILATERAL, 1, 4, 1,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 1, DM_POLYTOPE_QUADRILATERAL, 1, 1, 1, DM_POLYTOPE_QUADRILATERAL, 1, 2, 1, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1,  3, 0, DM_POLYTOPE_QUADRILATERAL, 0,    0,
                                     DM_POLYTOPE_QUADRILATERAL, 1, 0, 2, DM_POLYTOPE_QUADRILATERAL, 1, 1, 2, DM_POLYTOPE_QUADRILATERAL, 0,    2, DM_POLYTOPE_QUADRILATERAL, 1, 3, 1, DM_POLYTOPE_QUADRILATERAL, 0,    1, DM_POLYTOPE_QUADRILATERAL, 1, 4, 0};
  static PetscInt       ctripO[]  = {0, 0,
                                     0, 0, -2, -2,
                                     0, 0, -2, -2,
                                     0, 0, -2, -2,
                                    -4, 0, 0, -1,  0,  1,
                                    -4, 0, 0,  0,  0, -4,
                                    -4, 0, 0,  0, -1,  1};
  /* Add 1 edge and 4 quads inside every tensor quad prism, making 4 new hexahedra. */
  static DMPolytopeType tquadpT[]  = {DM_POLYTOPE_POINT_PRISM_TENSOR, DM_POLYTOPE_SEG_PRISM_TENSOR, DM_POLYTOPE_QUAD_PRISM_TENSOR};
  static PetscInt       tquadpS[]  = {1, 4, 4};
  static PetscInt       tquadpC[]  = {DM_POLYTOPE_POINT, 1, 0, 0, DM_POLYTOPE_POINT, 1, 1, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 0, DM_POLYTOPE_SEGMENT, 1, 1, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 1, DM_POLYTOPE_SEGMENT, 1, 1, 1, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 2, DM_POLYTOPE_SEGMENT, 1, 1, 2, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 4, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_SEGMENT, 1, 0, 3, DM_POLYTOPE_SEGMENT, 1, 1, 3, DM_POLYTOPE_POINT_PRISM_TENSOR, 1, 5, 0, DM_POLYTOPE_POINT_PRISM_TENSOR, 0, 0,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 0, DM_POLYTOPE_QUADRILATERAL, 1, 1, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    3, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 5, 1,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 1, DM_POLYTOPE_QUADRILATERAL, 1, 1, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 2, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    1, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    0,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 2, DM_POLYTOPE_QUADRILATERAL, 1, 1, 2, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 3, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 0, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2,
                                      DM_POLYTOPE_QUADRILATERAL, 1, 0, 3, DM_POLYTOPE_QUADRILATERAL, 1, 1, 3, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    3, DM_POLYTOPE_SEG_PRISM_TENSOR, 0,    2, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 4, 1, DM_POLYTOPE_SEG_PRISM_TENSOR, 1, 5, 0};
  static PetscInt       tquadpO[]  = {0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0, 0, 0,
                                      0, 0,  0,  0, -1,  0,
                                      0, 0,  0,  0,  0, -1,
                                      0, 0, -1,  0,  0,  0,
                                      0, 0,  0, -1,  0,  0};
  PetscBool convertTensor = PETSC_TRUE;

  PetscFunctionBeginHot;
  if (convertTensor) {
    switch (source) {
      case DM_POLYTOPE_POINT:
      case DM_POLYTOPE_SEGMENT:
      case DM_POLYTOPE_QUADRILATERAL:
      case DM_POLYTOPE_HEXAHEDRON:
        ierr = DMPlexCellRefinerRefine_Regular(cr, source, Nt, target, size, cone, ornt);CHKERRQ(ierr);
        break;
      case DM_POLYTOPE_POINT_PRISM_TENSOR: *Nt = 1; *target = tedgeT;  *size = tedgeS;  *cone = tedgeC;  *ornt = tedgeO;  break;
      case DM_POLYTOPE_SEG_PRISM_TENSOR:   *Nt = 2; *target = tquadT;  *size = tquadS;  *cone = tquadC;  *ornt = tquadO;  break;
      case DM_POLYTOPE_TRI_PRISM_TENSOR:   *Nt = 3; *target = ctripT;  *size = ctripS;  *cone = ctripC;  *ornt = ctripO;  break;
      case DM_POLYTOPE_QUAD_PRISM_TENSOR:  *Nt = 3; *target = tquadpT; *size = tquadpS; *cone = tquadpC; *ornt = tquadpO; break;
      case DM_POLYTOPE_TRIANGLE:           *Nt = 3; *target = triT;    *size = triS;    *cone = triC;    *ornt = triO;    break;
      case DM_POLYTOPE_TETRAHEDRON:        *Nt = 4; *target = tetT;    *size = tetS;    *cone = tetC;    *ornt = tetO;    break;
      case DM_POLYTOPE_TRI_PRISM:          *Nt = 4; *target = tripT;   *size = tripS;   *cone = tripC;   *ornt = tripO;   break;
      default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "No refinement strategy for %s", DMPolytopeTypes[source]);
    }
  } else {
    switch (source) {
      case DM_POLYTOPE_POINT:
      case DM_POLYTOPE_POINT_PRISM_TENSOR:
      case DM_POLYTOPE_SEGMENT:
      case DM_POLYTOPE_QUADRILATERAL:
      case DM_POLYTOPE_SEG_PRISM_TENSOR:
      case DM_POLYTOPE_HEXAHEDRON:
      case DM_POLYTOPE_QUAD_PRISM_TENSOR:
        ierr = DMPlexCellRefinerRefine_Regular(cr, source, Nt, target, size, cone, ornt);CHKERRQ(ierr);
        break;
      case DM_POLYTOPE_TRIANGLE:           *Nt = 3; *target = triT;    *size = triS;    *cone = triC;    *ornt = triO;    break;
      case DM_POLYTOPE_TETRAHEDRON:        *Nt = 4; *target = tetT;    *size = tetS;    *cone = tetC;    *ornt = tetO;    break;
      case DM_POLYTOPE_TRI_PRISM:          *Nt = 4; *target = tripT;   *size = tripS;   *cone = tripC;   *ornt = tripO;   break;
      case DM_POLYTOPE_TRI_PRISM_TENSOR:   *Nt = 3; *target = ttripT;  *size = ttripS;  *cone = ttripC;  *ornt = ttripO;  break;
      default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "No refinement strategy for %s", DMPolytopeTypes[source]);
    }
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerRefine_ToSimplex(DMPlexCellRefiner cr, DMPolytopeType source, PetscInt *Nt, DMPolytopeType *target[], PetscInt *size[], PetscInt *cone[], PetscInt *ornt[])
{
  PetscErrorCode ierr;

  PetscFunctionBeginHot;
  switch (source) {
    case DM_POLYTOPE_POINT:
    case DM_POLYTOPE_SEGMENT:
    case DM_POLYTOPE_POINT_PRISM_TENSOR:
    case DM_POLYTOPE_TRIANGLE:
    case DM_POLYTOPE_TETRAHEDRON:
    case DM_POLYTOPE_TRI_PRISM:
    case DM_POLYTOPE_TRI_PRISM_TENSOR:
    case DM_POLYTOPE_QUADRILATERAL:
    case DM_POLYTOPE_SEG_PRISM_TENSOR:
    case DM_POLYTOPE_HEXAHEDRON:
    case DM_POLYTOPE_QUAD_PRISM_TENSOR:
      ierr = DMPlexCellRefinerRefine_Regular(cr, source, Nt, target, size, cone, ornt);CHKERRQ(ierr);
      break;
    default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "No refinement strategy for %s", DMPolytopeTypes[source]);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode CellRefinerCreateOffset_Internal(DMPlexCellRefiner cr, PetscInt ctOrder[], PetscInt ctStart[], PetscInt **offset)
{
  PetscInt       c, cN, *off;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscCalloc1(DM_NUM_POLYTOPES*DM_NUM_POLYTOPES, &off);CHKERRQ(ierr);
  for (c = DM_POLYTOPE_POINT; c < DM_NUM_POLYTOPES; ++c) {
    const DMPolytopeType ct = (DMPolytopeType) c;
    for (cN = DM_POLYTOPE_POINT; cN < DM_NUM_POLYTOPES; ++cN) {
      const DMPolytopeType ctNew = (DMPolytopeType) cN;
      DMPolytopeType      *rct;
      PetscInt            *rsize, *cone, *ornt;
      PetscInt             Nct, n, i;

      if (DMPolytopeTypeGetDim(ct) < 0 || DMPolytopeTypeGetDim(ctNew) < 0) {off[ct*DM_NUM_POLYTOPES+ctNew] = -1; break;}
      off[ct*DM_NUM_POLYTOPES+ctNew] = 0;
      for (i = DM_POLYTOPE_POINT; i < DM_NUM_POLYTOPES; ++i) {
        const DMPolytopeType ict  = (DMPolytopeType) ctOrder[i];
        const DMPolytopeType ictn = (DMPolytopeType) ctOrder[i+1];

        ierr = DMPlexCellRefinerRefine(cr, ict, &Nct, &rct, &rsize, &cone, &ornt);CHKERRQ(ierr);
        if (ict == ct) {
          for (n = 0; n < Nct; ++n) if (rct[n] == ctNew) break;
          if (n == Nct) off[ct*DM_NUM_POLYTOPES+ctNew] = -1;
          break;
        }
        for (n = 0; n < Nct; ++n) if (rct[n] == ctNew) off[ct*DM_NUM_POLYTOPES+ctNew] += (ctStart[ictn]-ctStart[ict]) * rsize[n];
      }
    }
  }
  *offset = off;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerSetStarts(DMPlexCellRefiner cr, const PetscInt ctStart[], const PetscInt ctStartNew[])
{
  const PetscInt ctSize = DM_NUM_POLYTOPES+1;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (cr->setupcalled) SETERRQ(PetscObjectComm((PetscObject) cr), PETSC_ERR_ARG_WRONGSTATE, "Must call this function before DMPlexCellRefinerSetUp()");
  ierr = PetscCalloc2(ctSize, &cr->ctStart, ctSize, &cr->ctStartNew);CHKERRQ(ierr);
  ierr = PetscArraycpy(cr->ctStart,    ctStart,    ctSize);CHKERRQ(ierr);
  ierr = PetscArraycpy(cr->ctStartNew, ctStartNew, ctSize);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Construct cell type order since we must loop over cell types in depth order */
PetscErrorCode DMPlexCreateCellTypeOrder_Internal(DMPolytopeType ctCell, PetscInt *ctOrder[], PetscInt *ctOrderInv[])
{
  PetscInt      *ctO, *ctOInv;
  PetscInt       c, d, off = 0;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscCalloc2(DM_NUM_POLYTOPES+1, &ctO, DM_NUM_POLYTOPES+1, &ctOInv);CHKERRQ(ierr);
  for (d = 3; d >= DMPolytopeTypeGetDim(ctCell); --d) {
    for (c = 0; c <= DM_NUM_POLYTOPES; ++c) {
      if (DMPolytopeTypeGetDim((DMPolytopeType) c) != d) continue;
      ctO[off++] = c;
    }
  }
  if (DMPolytopeTypeGetDim(ctCell) != 0) {
    for (c = 0; c <= DM_NUM_POLYTOPES; ++c) {
      if (DMPolytopeTypeGetDim((DMPolytopeType) c) != 0) continue;
      ctO[off++] = c;
    }
  }
  for (d = DMPolytopeTypeGetDim(ctCell)-1; d > 0; --d) {
    for (c = 0; c <= DM_NUM_POLYTOPES; ++c) {
      if (DMPolytopeTypeGetDim((DMPolytopeType) c) != d) continue;
      ctO[off++] = c;
    }
  }
  for (c = 0; c <= DM_NUM_POLYTOPES; ++c) {
    if (DMPolytopeTypeGetDim((DMPolytopeType) c) >= 0) continue;
    ctO[off++] = c;
  }
  if (off != DM_NUM_POLYTOPES+1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Invalid offset %D for cell type order", off);
  for (c = 0; c <= DM_NUM_POLYTOPES; ++c) {
    ctOInv[ctO[c]] = c;
  }
  *ctOrder    = ctO;
  *ctOrderInv = ctOInv;
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlexCellRefinerSetUp(DMPlexCellRefiner cr)
{
  DM             dm = cr->dm;
  DMPolytopeType ctCell;
  PetscInt       pStart, pEnd, p, c;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeader(cr, 1);
  if (cr->setupcalled) PetscFunctionReturn(0);
  ierr = DMPlexGetChart(dm, &pStart, &pEnd);CHKERRQ(ierr);
  if (pEnd > pStart) {ierr = DMPlexGetCellType(dm, 0, &ctCell);CHKERRQ(ierr);}
  else               {
    PetscInt dim;
    ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
    switch (dim) {
      case 0: ctCell = DM_POLYTOPE_POINT;break;
      case 1: ctCell = DM_POLYTOPE_SEGMENT;break;
      case 2: ctCell = DM_POLYTOPE_TRIANGLE;break;
      case 3: ctCell = DM_POLYTOPE_TETRAHEDRON;break;
      default: ctCell = DM_POLYTOPE_TETRAHEDRON;
    }
  }
  ierr = DMPlexCreateCellTypeOrder_Internal(ctCell, &cr->ctOrder, &cr->ctOrderInv);CHKERRQ(ierr);
  /* Construct sizes and offsets for each cell type */
  if (!cr->ctStart) {
    PetscInt *ctS, *ctSN, *ctC, *ctCN;

    ierr = PetscCalloc2(DM_NUM_POLYTOPES+1, &ctS, DM_NUM_POLYTOPES+1, &ctSN);CHKERRQ(ierr);
    ierr = PetscCalloc2(DM_NUM_POLYTOPES+1, &ctC, DM_NUM_POLYTOPES+1, &ctCN);CHKERRQ(ierr);
    for (p = pStart; p < pEnd; ++p) {
      DMPolytopeType  ct;
      DMPolytopeType *rct;
      PetscInt       *rsize, *cone, *ornt;
      PetscInt        Nct, n;

      ierr = DMPlexGetCellType(dm, p, &ct);CHKERRQ(ierr);
      if ((PetscInt) ct < 0) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "No cell type for point %D", p);
      ++ctC[ct];
      ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &cone, &ornt);CHKERRQ(ierr);
      for (n = 0; n < Nct; ++n) ctCN[rct[n]] += rsize[n];
    }
    for (c = 0; c < DM_NUM_POLYTOPES; ++c) {
      const PetscInt ct  = cr->ctOrder[c];
      const PetscInt ctn = cr->ctOrder[c+1];

      ctS[ctn]  = ctS[ct]  + ctC[ct];
      ctSN[ctn] = ctSN[ct] + ctCN[ct];
    }
    ierr = PetscFree2(ctC, ctCN);CHKERRQ(ierr);
    cr->ctStart    = ctS;
    cr->ctStartNew = ctSN;
  }
  ierr = CellRefinerCreateOffset_Internal(cr, cr->ctOrder, cr->ctStart, &cr->offset);CHKERRQ(ierr);
  cr->setupcalled = PETSC_TRUE;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerView_Ascii(DMPlexCellRefiner cr, PetscViewer v)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscViewerASCIIPrintf(v, "Cell Refiner: %s\n", DMPlexCellRefinerTypes[cr->type]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
  DMPlexCellRefinerView - Views a DMPlexCellRefiner object

  Collective on cr

  Input Parameters:
+ cr     - The DMPlexCellRefiner object
- viewer - The PetscViewer object

  Level: beginner

.seealso: DMPlexCellRefinerCreate()
*/
static PetscErrorCode DMPlexCellRefinerView(DMPlexCellRefiner cr, PetscViewer viewer)
{
  PetscBool      iascii;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeader(cr, 1);
  if (viewer) PetscValidHeaderSpecific(viewer, PETSC_VIEWER_CLASSID, 2);
  if (!viewer) {ierr = PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject) cr), &viewer);CHKERRQ(ierr);}
  ierr = PetscObjectTypeCompare((PetscObject) viewer, PETSCVIEWERASCII, &iascii);CHKERRQ(ierr);
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  if (iascii) {ierr = DMPlexCellRefinerView_Ascii(cr, viewer);CHKERRQ(ierr);}
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlexCellRefinerDestroy(DMPlexCellRefiner *cr)
{
  PetscInt       c;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!*cr) PetscFunctionReturn(0);
  PetscValidHeaderSpecific(*cr, DM_CLASSID, 1);
  ierr = PetscObjectDereference((PetscObject) (*cr)->dm);CHKERRQ(ierr);
  ierr = PetscFree2((*cr)->ctOrder, (*cr)->ctOrderInv);CHKERRQ(ierr);
  ierr = PetscFree2((*cr)->ctStart, (*cr)->ctStartNew);CHKERRQ(ierr);
  ierr = PetscFree((*cr)->offset);CHKERRQ(ierr);
  for (c = 0; c < DM_NUM_POLYTOPES; ++c) {
    ierr = PetscFEDestroy(&(*cr)->coordFE[c]);CHKERRQ(ierr);
    ierr = PetscFEGeomDestroy(&(*cr)->refGeom[c]);CHKERRQ(ierr);
  }
  ierr = PetscFree2((*cr)->coordFE, (*cr)->refGeom);CHKERRQ(ierr);
  ierr = PetscHeaderDestroy(cr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMPlexCellRefinerCreate(DM dm, DMPlexCellRefiner *cr)
{
  DMPlexCellRefiner tmp;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  PetscValidPointer(cr, 2);
  *cr  = NULL;
  ierr = PetscHeaderCreate(tmp, DM_CLASSID, "DMPlexCellRefiner", "Cell Refiner", "DMPlexCellRefiner", PETSC_COMM_SELF, DMPlexCellRefinerDestroy, DMPlexCellRefinerView);CHKERRQ(ierr);
  tmp->setupcalled = PETSC_FALSE;

  tmp->dm = dm;
  ierr = PetscObjectReference((PetscObject) dm);CHKERRQ(ierr);
  ierr = DMPlexGetCellRefinerType(dm, &tmp->type);CHKERRQ(ierr);
  switch (tmp->type) {
    case DM_REFINER_REGULAR:
      tmp->ops->refine                  = DMPlexCellRefinerRefine_Regular;
      tmp->ops->mapsubcells             = DMPlexCellRefinerMapSubcells_Regular;
      tmp->ops->getcellvertices         = DMPlexCellRefinerGetCellVertices_Regular;
      tmp->ops->getsubcellvertices      = DMPlexCellRefinerGetSubcellVertices_Regular;
      tmp->ops->getaffinetransforms     = DMPlexCellRefinerGetAffineTransforms_Regular;
      tmp->ops->getaffinefacetransforms = DMPlexCellRefinerGetAffineFaceTransforms_Regular;
      break;
    case DM_REFINER_TO_BOX:
      tmp->ops->refine             = DMPlexCellRefinerRefine_ToBox;
      tmp->ops->mapsubcells        = DMPlexCellRefinerMapSubcells_ToBox;
      tmp->ops->getcellvertices    = DMPlexCellRefinerGetCellVertices_ToBox;
      tmp->ops->getsubcellvertices = DMPlexCellRefinerGetSubcellVertices_ToBox;
      break;
    case DM_REFINER_TO_SIMPLEX:
      tmp->ops->refine      = DMPlexCellRefinerRefine_ToSimplex;
      tmp->ops->mapsubcells = DMPlexCellRefinerMapSubcells_ToSimplex;
      break;
    default: SETERRQ1(PetscObjectComm((PetscObject) dm), PETSC_ERR_ARG_WRONG, "Invalid cell refiner type %s", DMPlexCellRefinerTypes[tmp->type]);
  }
  ierr = PetscCalloc2(DM_NUM_POLYTOPES, &tmp->coordFE, DM_NUM_POLYTOPES, &tmp->refGeom);CHKERRQ(ierr);
  *cr = tmp;
  PetscFunctionReturn(0);
}

/*@
  DMPlexCellRefinerGetAffineTransforms - Gets the affine map from the reference cell to each subcell

  Input Parameters:
+ cr - The DMPlexCellRefiner object
- ct - The cell type

  Output Parameters:
+ Nc   - The number of subcells produced from this cell type
. v0   - The translation of the first vertex for each subcell
. J    - The Jacobian for each subcell (map from reference cell to subcell)
- invJ - The inverse Jacobian for each subcell

  Level: developer

.seealso: DMPlexCellRefinerGetAffineFaceTransforms(), Create()
@*/
PetscErrorCode DMPlexCellRefinerGetAffineTransforms(DMPlexCellRefiner cr, DMPolytopeType ct, PetscInt *Nc, PetscReal *v0[], PetscReal *J[], PetscReal *invJ[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!cr->ops->getaffinetransforms) SETERRQ(PetscObjectComm((PetscObject) cr), PETSC_ERR_SUP, "No support for affine transforms from this refiner");
  ierr = (*cr->ops->getaffinetransforms)(cr, ct, Nc, v0, J, invJ);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
  DMPlexCellRefinerGetAffineFaceTransforms - Gets the affine map from the reference face cell to each face in the given cell

  Input Parameters:
+ cr - The DMPlexCellRefiner object
- ct - The cell type

  Output Parameters:
+ Nf   - The number of faces for this cell type
. v0   - The translation of the first vertex for each face
. J    - The Jacobian for each face (map from original cell to subcell)
. invJ - The inverse Jacobian for each face
- detJ - The determinant of the Jacobian for each face

  Note: The Jacobian and inverse Jacboian will be rectangular, and the inverse is really a generalized inverse.

  Level: developer

.seealso: DMPlexCellRefinerGetAffineTransforms(), Create()
@*/
PetscErrorCode DMPlexCellRefinerGetAffineFaceTransforms(DMPlexCellRefiner cr, DMPolytopeType ct, PetscInt *Nf, PetscReal *v0[], PetscReal *J[], PetscReal *invJ[], PetscReal *detJ[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!cr->ops->getaffinefacetransforms) SETERRQ(PetscObjectComm((PetscObject) cr), PETSC_ERR_SUP, "No support for affine face transforms from this refiner");
  ierr = (*cr->ops->getaffinefacetransforms)(cr, ct, Nf, v0, J, invJ, detJ);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* Numbering regularly refined meshes

   We want the numbering of the new mesh to respect the same depth stratification as the old mesh. We first compute
   the number of new points at each depth. This means that offsets for each depth can be computed, making no assumptions
   about the order of different cell types.

   However, when we want to order different depth strata, it will be very useful to make assumptions about contiguous
   numbering of different cell types, especially if we want to compute new numberings without communication. Therefore, we
   will require that cells are numbering contiguously for each cell type, and that those blocks come in the same order as
   the cell type enumeration within a given depth stratum.

   Thus, at each depth, each cell type will add a certain number of points at that depth. To get the new point number, we
   start at the new depth offset, run through all prior cell types incrementing by the total addition from that type, then
   offset by the old cell type number and replica number for the insertion.
*/
PetscErrorCode DMPlexCellRefinerGetNewPoint(DMPlexCellRefiner cr, DMPolytopeType ct, DMPolytopeType ctNew, PetscInt p, PetscInt r, PetscInt *pNew)
{
  DMPolytopeType  *rct;
  PetscInt        *rsize, *cone, *ornt;
  PetscInt         Nct, n;
  PetscInt         off  = cr->offset[ct*DM_NUM_POLYTOPES+ctNew];
  PetscInt         ctS  = cr->ctStart[ct],       ctE  = cr->ctStart[cr->ctOrder[cr->ctOrderInv[ct]+1]];
  PetscInt         ctSN = cr->ctStartNew[ctNew], ctEN = cr->ctStartNew[cr->ctOrder[cr->ctOrderInv[ctNew]+1]];
  PetscInt         newp = ctSN;
  PetscErrorCode   ierr;

  PetscFunctionBeginHot;
  if ((p < ctS) || (p >= ctE)) SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Point %D is not a %s [%D, %D)", p, DMPolytopeTypes[ct], ctS, ctE);
  if (off < 0) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Cell type %s does not produce type %s", DMPolytopeTypes[ct], DMPolytopeTypes[ctNew]);

  newp += off;
  ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &cone, &ornt);CHKERRQ(ierr);
  for (n = 0; n < Nct; ++n) {
    if (rct[n] == ctNew) {
      if (rsize[n] && r >= rsize[n])
        SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Replica number %D should be in [0, %D) for subcell type %s in cell type %s", r, rsize[n], DMPolytopeTypes[rct[n]], DMPolytopeTypes[ct]);
      newp += (p - ctS) * rsize[n] + r;
      break;
    }
  }

  if ((newp < ctSN) || (newp >= ctEN)) SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_PLIB, "New point %D is not a %s [%D, %D)", newp, DMPolytopeTypes[ctNew], ctSN, ctEN);
  *pNew = newp;
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerSetConeSizes(DMPlexCellRefiner cr, DM rdm)
{
  DM              dm = cr->dm;
  PetscInt        pStart, pEnd, p, pNew;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  /* Must create the celltype label here so that we do not automatically try to compute the types */
  ierr = DMCreateLabel(rdm, "celltype");CHKERRQ(ierr);
  ierr = DMPlexGetChart(dm, &pStart, &pEnd);CHKERRQ(ierr);
  for (p = pStart; p < pEnd; ++p) {
    DMPolytopeType  ct;
    DMPolytopeType *rct;
    PetscInt       *rsize, *rcone, *rornt;
    PetscInt        Nct, n, r;

    ierr = DMPlexGetCellType(dm, p, &ct);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
    for (n = 0; n < Nct; ++n) {
      for (r = 0; r < rsize[n]; ++r) {
        ierr = DMPlexCellRefinerGetNewPoint(cr, ct, rct[n], p, r, &pNew);CHKERRQ(ierr);
        ierr = DMPlexSetConeSize(rdm, pNew, DMPolytopeTypeGetConeSize(rct[n]));CHKERRQ(ierr);
        ierr = DMPlexSetCellType(rdm, pNew, rct[n]);CHKERRQ(ierr);
      }
    }
  }
  {
    DMLabel  ctLabel;
    DM_Plex *plex = (DM_Plex *) rdm->data;

    ierr = DMPlexGetCellTypeLabel(rdm, &ctLabel);CHKERRQ(ierr);
    ierr = PetscObjectStateGet((PetscObject) ctLabel, &plex->celltypeState);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerSetCones(DMPlexCellRefiner cr, DM rdm)
{
  DM             dm = cr->dm;
  DMPolytopeType ct;
  PetscInt      *coneNew, *orntNew;
  PetscInt       maxConeSize = 0, pStart, pEnd, p, pNew;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  for (p = 0; p < DM_NUM_POLYTOPES; ++p) maxConeSize = PetscMax(maxConeSize, DMPolytopeTypeGetConeSize((DMPolytopeType) p));
  ierr = PetscMalloc2(maxConeSize, &coneNew, maxConeSize, &orntNew);CHKERRQ(ierr);
  ierr = DMPlexGetChart(dm, &pStart, &pEnd);CHKERRQ(ierr);
  for (p = pStart; p < pEnd; ++p) {
    const PetscInt *cone, *ornt;
    PetscInt        coff, ooff, c;
    DMPolytopeType *rct;
    PetscInt       *rsize, *rcone, *rornt;
    PetscInt        Nct, n, r;

    ierr = DMPlexGetCellType(dm, p, &ct);CHKERRQ(ierr);
    ierr = DMPlexGetCone(dm, p, &cone);CHKERRQ(ierr);
    ierr = DMPlexGetConeOrientation(dm, p, &ornt);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
    for (n = 0, coff = 0, ooff = 0; n < Nct; ++n) {
      const DMPolytopeType ctNew    = rct[n];
      const PetscInt       csizeNew = DMPolytopeTypeGetConeSize(ctNew);

      for (r = 0; r < rsize[n]; ++r) {
        /* pNew is a subcell produced by subdividing p */
        ierr = DMPlexCellRefinerGetNewPoint(cr, ct, rct[n], p, r, &pNew);CHKERRQ(ierr);
        for (c = 0; c < csizeNew; ++c) {
          PetscInt             ppp   = -1;                             /* Parent Parent point: Parent of point pp */
          PetscInt             pp    = p;                              /* Parent point: Point in the original mesh producing new cone point */
          PetscInt             po    = 0;                              /* Orientation of parent point pp in parent parent point ppp */
          DMPolytopeType       pct   = ct;                             /* Parent type:  Cell type for parent of new cone point */
          const PetscInt      *pcone = cone;                           /* Parent cone: Cone of parent point pp */
          PetscInt             pr    = -1;                             /* Replica number of pp that produces new cone point  */
          const DMPolytopeType ft    = (DMPolytopeType) rcone[coff++]; /* Cell type for new cone point of pNew */
          const PetscInt       fn    = rcone[coff++];                  /* Number of cones of p that need to be taken when producing new cone point */
          PetscInt             fo    = rornt[ooff++];                  /* Orientation of new cone point in pNew */
          PetscInt             lc;

          /* Get the type (pct) and point number (pp) of the parent point in the original mesh which produces this cone point */
          for (lc = 0; lc < fn; ++lc) {
            const PetscInt *ppornt;
            PetscInt        pcp;

            ierr = DMPolytopeMapCell(pct, po, rcone[coff++], &pcp);CHKERRQ(ierr);
            ppp  = pp;
            pp   = pcone[pcp];
            ierr = DMPlexGetCellType(dm, pp, &pct);CHKERRQ(ierr);
            ierr = DMPlexGetCone(dm, pp, &pcone);CHKERRQ(ierr);
            ierr = DMPlexGetConeOrientation(dm, ppp, &ppornt);CHKERRQ(ierr);
            po   = ppornt[pcp];
          }
          pr = rcone[coff++];
          /* Orientation po of pp maps (pr, fo) -> (pr', fo') */
          ierr = DMPlexCellRefinerMapSubcells(cr, pct, po, ft, pr, fo, &pr, &fo);CHKERRQ(ierr);
          ierr = DMPlexCellRefinerGetNewPoint(cr, pct, ft, pp, pr, &coneNew[c]);CHKERRQ(ierr);
          orntNew[c] = fo;
        }
        ierr = DMPlexSetCone(rdm, pNew, coneNew);CHKERRQ(ierr);
        ierr = DMPlexSetConeOrientation(rdm, pNew, orntNew);CHKERRQ(ierr);
      }
    }
  }
  ierr = PetscFree2(coneNew, orntNew);CHKERRQ(ierr);
  ierr = DMPlexSymmetrize(rdm);CHKERRQ(ierr);
  ierr = DMPlexStratify(rdm);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerGetCoordinateFE(DMPlexCellRefiner cr, DMPolytopeType ct, PetscFE *fe)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!cr->coordFE[ct]) {
    PetscInt  dim, cdim;
    PetscBool isSimplex;

    switch (ct) {
      case DM_POLYTOPE_SEGMENT:       dim = 1; isSimplex = PETSC_TRUE;  break;
      case DM_POLYTOPE_TRIANGLE:      dim = 2; isSimplex = PETSC_TRUE;  break;
      case DM_POLYTOPE_QUADRILATERAL: dim = 2; isSimplex = PETSC_FALSE; break;
      case DM_POLYTOPE_TETRAHEDRON:   dim = 3; isSimplex = PETSC_TRUE;  break;
      case DM_POLYTOPE_HEXAHEDRON:    dim = 3; isSimplex = PETSC_FALSE; break;
      default: SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_SUP, "No coordinate FE for cell type %s", DMPolytopeTypes[ct]);
    }
    ierr = DMGetCoordinateDim(cr->dm, &cdim);CHKERRQ(ierr);
    ierr = PetscFECreateLagrange(PETSC_COMM_SELF, dim, cdim, isSimplex, 1, PETSC_DETERMINE, &cr->coordFE[ct]);CHKERRQ(ierr);
    {
      PetscDualSpace  dsp;
      PetscQuadrature quad;
      DM              K;
      PetscFEGeom    *cg;
      PetscReal      *Xq, *xq, *wq;
      PetscInt        Nq, q;

      ierr = DMPlexCellRefinerGetCellVertices(cr, ct, &Nq, &Xq);CHKERRQ(ierr);
      ierr = PetscMalloc1(Nq*cdim, &xq);CHKERRQ(ierr);
      for (q = 0; q < Nq*cdim; ++q) xq[q] = Xq[q];
      ierr = PetscMalloc1(Nq, &wq);CHKERRQ(ierr);
      for (q = 0; q < Nq; ++q) wq[q] = 1.0;
      ierr = PetscQuadratureCreate(PETSC_COMM_SELF, &quad);CHKERRQ(ierr);
      ierr = PetscQuadratureSetData(quad, dim, 1, Nq, xq, wq);CHKERRQ(ierr);
      ierr = PetscFESetQuadrature(cr->coordFE[ct], quad);CHKERRQ(ierr);

      ierr = PetscFEGetDualSpace(cr->coordFE[ct], &dsp);CHKERRQ(ierr);
      ierr = PetscDualSpaceGetDM(dsp, &K);CHKERRQ(ierr);
      ierr = PetscFEGeomCreate(quad, 1, cdim, PETSC_FALSE, &cr->refGeom[ct]);CHKERRQ(ierr);
      cg   = cr->refGeom[ct];
      ierr = DMPlexComputeCellGeometryFEM(K, 0, NULL, cg->v, cg->J, cg->invJ, cg->detJ);CHKERRQ(ierr);
      ierr = PetscQuadratureDestroy(&quad);CHKERRQ(ierr);
    }
  }
  *fe = cr->coordFE[ct];
  PetscFunctionReturn(0);
}

/*
  DMPlexCellRefinerMapLocalizedCoordinates - Given a cell of type ct with localized coordinates x, we generate localized coordinates xr for subcell r of type rct.

  Not collective

  Input Parameters:
+ cr  - The DMPlexCellRefiner
. ct  - The type of the parent cell
. rct - The type of the produced cell
. r   - The index of the produced cell
- x   - The localized coordinates for the parent cell

  Output Parameter:
. xr  - The localized coordinates for the produced cell

  Level: developer

.seealso: DMPlexCellRefinerSetCoordinates()
*/
static PetscErrorCode DMPlexCellRefinerMapLocalizedCoordinates(DMPlexCellRefiner cr, DMPolytopeType ct, DMPolytopeType rct, PetscInt r, const PetscScalar x[], PetscScalar xr[])
{
  PetscFE        fe = NULL;
  PetscInt       cdim, Nv, v, *subcellV;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexCellRefinerGetCoordinateFE(cr, ct, &fe);CHKERRQ(ierr);
  ierr = DMPlexCellRefinerGetSubcellVertices(cr, ct, rct, r, &Nv, &subcellV);CHKERRQ(ierr);
  ierr = PetscFEGetNumComponents(fe, &cdim);CHKERRQ(ierr);
  for (v = 0; v < Nv; ++v) {
    ierr = PetscFEInterpolate_Static(fe, x, cr->refGeom[ct], subcellV[v], &xr[v*cdim]);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerSetCoordinates(DMPlexCellRefiner cr, DM rdm)
{
  DM                    dm = cr->dm, cdm;
  PetscSection          coordSection, coordSectionNew;
  Vec                   coordsLocal, coordsLocalNew;
  const PetscScalar    *coords;
  PetscScalar          *coordsNew;
  const DMBoundaryType *bd;
  const PetscReal      *maxCell, *L;
  PetscBool             isperiodic, localizeVertices = PETSC_FALSE, localizeCells = PETSC_FALSE;
  PetscInt              dE, d, cStart, cEnd, c, vStartNew, vEndNew, v, pStart, pEnd, p, ocStart, ocEnd;
  PetscErrorCode        ierr;

  PetscFunctionBegin;
  ierr = DMGetCoordinateDM(dm, &cdm);CHKERRQ(ierr);
  ierr = DMGetPeriodicity(dm, &isperiodic, &maxCell, &L, &bd);CHKERRQ(ierr);
  /* Determine if we need to localize coordinates when generating them */
  if (isperiodic) {
    localizeVertices = PETSC_TRUE;
    if (!maxCell) {
      PetscBool localized;
      ierr = DMGetCoordinatesLocalized(dm, &localized);CHKERRQ(ierr);
      if (!localized) SETERRQ(PetscObjectComm((PetscObject) dm), PETSC_ERR_USER, "Cannot refine a periodic mesh if coordinates have not been localized");
      localizeCells = PETSC_TRUE;
    }
  }

  ierr = DMGetCoordinateSection(dm, &coordSection);CHKERRQ(ierr);
  ierr = PetscSectionGetFieldComponents(coordSection, 0, &dE);CHKERRQ(ierr);
  if (maxCell) {
    PetscReal maxCellNew[3];

    for (d = 0; d < dE; ++d) maxCellNew[d] = maxCell[d]/2.0;
    ierr = DMSetPeriodicity(rdm, isperiodic, maxCellNew, L, bd);CHKERRQ(ierr);
  } else {
    ierr = DMSetPeriodicity(rdm, isperiodic, maxCell, L, bd);CHKERRQ(ierr);
  }
  ierr = PetscSectionCreate(PetscObjectComm((PetscObject) dm), &coordSectionNew);CHKERRQ(ierr);
  ierr = PetscSectionSetNumFields(coordSectionNew, 1);CHKERRQ(ierr);
  ierr = PetscSectionSetFieldComponents(coordSectionNew, 0, dE);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(rdm, 0, &vStartNew, &vEndNew);CHKERRQ(ierr);
  if (localizeCells) {ierr = PetscSectionSetChart(coordSectionNew, 0,         vEndNew);CHKERRQ(ierr);}
  else               {ierr = PetscSectionSetChart(coordSectionNew, vStartNew, vEndNew);CHKERRQ(ierr);}

  /* Localization should be inherited */
  /*   Stefano calculates parent cells for each new cell for localization */
  /*   Localized cells need coordinates of closure */
  for (v = vStartNew; v < vEndNew; ++v) {
    ierr = PetscSectionSetDof(coordSectionNew, v, dE);CHKERRQ(ierr);
    ierr = PetscSectionSetFieldDof(coordSectionNew, v, 0, dE);CHKERRQ(ierr);
  }
  if (localizeCells) {
    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
    for (c = cStart; c < cEnd; ++c) {
      PetscInt dof;

      ierr = PetscSectionGetDof(coordSection, c, &dof); CHKERRQ(ierr);
      if (dof) {
        DMPolytopeType  ct;
        DMPolytopeType *rct;
        PetscInt       *rsize, *rcone, *rornt;
        PetscInt        dim, cNew, Nct, n, r;

        ierr = DMPlexGetCellType(dm, c, &ct);CHKERRQ(ierr);
        dim  = DMPolytopeTypeGetDim(ct);
        ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
        /* This allows for different cell types */
        for (n = 0; n < Nct; ++n) {
          if (dim != DMPolytopeTypeGetDim(rct[n])) continue;
          for (r = 0; r < rsize[n]; ++r) {
            PetscInt *closure = NULL;
            PetscInt  clSize, cl, Nv = 0;

            ierr = DMPlexCellRefinerGetNewPoint(cr, ct, rct[n], c, r, &cNew);CHKERRQ(ierr);
            ierr = DMPlexGetTransitiveClosure(rdm, cNew, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
            for (cl = 0; cl < clSize*2; cl += 2) {if ((closure[cl] >= vStartNew) && (closure[cl] < vEndNew)) ++Nv;}
            ierr = DMPlexRestoreTransitiveClosure(rdm, cNew, PETSC_TRUE, &clSize, &closure);CHKERRQ(ierr);
            ierr = PetscSectionSetDof(coordSectionNew, cNew, Nv * dE);CHKERRQ(ierr);
            ierr = PetscSectionSetFieldDof(coordSectionNew, cNew, 0, Nv * dE);CHKERRQ(ierr);
          }
        }
      }
    }
  }
  ierr = PetscSectionSetUp(coordSectionNew);CHKERRQ(ierr);
  ierr = DMViewFromOptions(dm, NULL, "-coarse_dm_view");CHKERRQ(ierr);
  ierr = DMSetCoordinateSection(rdm, PETSC_DETERMINE, coordSectionNew);CHKERRQ(ierr);
  {
    VecType     vtype;
    PetscInt    coordSizeNew, bs;
    const char *name;

    ierr = DMGetCoordinatesLocal(dm, &coordsLocal);CHKERRQ(ierr);
    ierr = VecCreate(PETSC_COMM_SELF, &coordsLocalNew);CHKERRQ(ierr);
    ierr = PetscSectionGetStorageSize(coordSectionNew, &coordSizeNew);CHKERRQ(ierr);
    ierr = VecSetSizes(coordsLocalNew, coordSizeNew, PETSC_DETERMINE);CHKERRQ(ierr);
    ierr = PetscObjectGetName((PetscObject) coordsLocal, &name);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) coordsLocalNew, name);CHKERRQ(ierr);
    ierr = VecGetBlockSize(coordsLocal, &bs);CHKERRQ(ierr);
    ierr = VecSetBlockSize(coordsLocalNew, bs);CHKERRQ(ierr);
    ierr = VecGetType(coordsLocal, &vtype);CHKERRQ(ierr);
    ierr = VecSetType(coordsLocalNew, vtype);CHKERRQ(ierr);
  }
  ierr = VecGetArrayRead(coordsLocal, &coords);CHKERRQ(ierr);
  ierr = VecGetArray(coordsLocalNew, &coordsNew);CHKERRQ(ierr);
  ierr = PetscSectionGetChart(coordSection, &ocStart, &ocEnd);CHKERRQ(ierr);
  ierr = DMPlexGetChart(dm, &pStart, &pEnd);CHKERRQ(ierr);
  /* First set coordinates for vertices*/
  for (p = pStart; p < pEnd; ++p) {
    DMPolytopeType  ct;
    DMPolytopeType *rct;
    PetscInt       *rsize, *rcone, *rornt;
    PetscInt        Nct, n, r;
    PetscBool       hasVertex = PETSC_FALSE, isLocalized = PETSC_FALSE;

    ierr = DMPlexGetCellType(dm, p, &ct);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
    for (n = 0; n < Nct; ++n) {
      if (rct[n] == DM_POLYTOPE_POINT) {hasVertex = PETSC_TRUE; break;}
    }
    if (localizeVertices && ct != DM_POLYTOPE_POINT && (p >= ocStart) && (p < ocEnd)) {
      PetscInt dof;
      ierr = PetscSectionGetDof(coordSection, p, &dof); CHKERRQ(ierr);
      if (dof) isLocalized = PETSC_TRUE;
    }
    if (hasVertex) {
      PetscScalar *pcoords = NULL;
      PetscScalar  vcoords[3] = {0., 0., 0.};
      PetscInt     Nc, Nv, v, d;

      ierr = DMPlexVecGetClosure(dm, coordSection, coordsLocal, p, &Nc, &pcoords);CHKERRQ(ierr);
      if (ct == DM_POLYTOPE_POINT) {
        for (d = 0; d < dE; ++d) vcoords[d] = pcoords[d];
      } else {
        if (localizeVertices) {
          PetscScalar anchor[3];

          for (d = 0; d < dE; ++d) anchor[d] = pcoords[d];
          if (!isLocalized) {
            Nv = Nc/dE;
            for (v = 0; v < Nv; ++v) {ierr = DMLocalizeAddCoordinate_Internal(dm, dE, anchor, &pcoords[v*dE], vcoords);CHKERRQ(ierr);}
          } else {
            Nv = Nc/(2*dE);
            for (v = Nv; v < Nv*2; ++v) {ierr = DMLocalizeAddCoordinate_Internal(dm, dE, anchor, &pcoords[v*dE], vcoords);CHKERRQ(ierr);}
          }
        } else {
          Nv = Nc/dE;
          for (v = 0; v < Nv; ++v) for (d = 0; d < dE; ++d) vcoords[d] += pcoords[v*dE+d];
        }
        for (d = 0; d < dE; ++d) vcoords[d] /= Nv;
      }
      ierr = DMPlexVecRestoreClosure(dm, coordSection, coordsLocal, p, &Nc, &pcoords);CHKERRQ(ierr);
      for (n = 0; n < Nct; ++n) {
        if (rct[n] != DM_POLYTOPE_POINT) continue;
        for (r = 0; r < rsize[n]; ++r) {
          PetscInt vNew, off;

          ierr = DMPlexCellRefinerGetNewPoint(cr, ct, rct[n], p, r, &vNew);CHKERRQ(ierr);
          ierr = PetscSectionGetOffset(coordSectionNew, vNew, &off);CHKERRQ(ierr);
          for (d = 0; d < dE; ++d) coordsNew[off+d] = vcoords[d];
        }
      }
    }
  }
  /* Then set coordinates for cells by localizing */
  for (p = pStart; p < pEnd; ++p) {
    DMPolytopeType  ct;
    DMPolytopeType *rct;
    PetscInt       *rsize, *rcone, *rornt;
    PetscInt        Nct, n, r;
    PetscBool       isLocalized = PETSC_FALSE;

    ierr = DMPlexGetCellType(dm, p, &ct);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
    if (localizeCells && ct != DM_POLYTOPE_POINT && (p >= ocStart) && (p < ocEnd)) {
      PetscInt dof;
      ierr = PetscSectionGetDof(coordSection, p, &dof); CHKERRQ(ierr);
      if (dof) isLocalized = PETSC_TRUE;
    }
    if (isLocalized) {
      const PetscScalar *pcoords;

      ierr = DMPlexPointLocalRead(cdm, p, coords, &pcoords);CHKERRQ(ierr);
      for (n = 0; n < Nct; ++n) {
        const PetscInt Nr = rsize[n];

        if (DMPolytopeTypeGetDim(ct) != DMPolytopeTypeGetDim(rct[n])) continue;
        for (r = 0; r < Nr; ++r) {
          PetscInt pNew, offNew;

          /* It looks like Stefano and Lisandro are allowing localized coordinates without defining the periodic boundary, which means that
             DMLocalizeCoordinate_Internal() will not work. Localized coordinates will have to have obtained by the affine map of the larger
             cell to the ones it produces. */
          ierr = DMPlexCellRefinerGetNewPoint(cr, ct, rct[n], p, r, &pNew);CHKERRQ(ierr);
          ierr = PetscSectionGetOffset(coordSectionNew, pNew, &offNew);CHKERRQ(ierr);
          ierr = DMPlexCellRefinerMapLocalizedCoordinates(cr, ct, rct[n], r, pcoords, &coordsNew[offNew]);CHKERRQ(ierr);
        }
      }
    }
  }
  ierr = VecRestoreArrayRead(coordsLocal, &coords);CHKERRQ(ierr);
  ierr = VecRestoreArray(coordsLocalNew, &coordsNew);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(rdm, coordsLocalNew);CHKERRQ(ierr);
  /* TODO Stefano has a final reduction if some hybrid coordinates cannot be found. (needcoords) Should not be needed. */
  ierr = VecDestroy(&coordsLocalNew);CHKERRQ(ierr);
  ierr = PetscSectionDestroy(&coordSectionNew);CHKERRQ(ierr);
  if (!localizeCells) {ierr = DMLocalizeCoordinates(rdm);CHKERRQ(ierr);}
  PetscFunctionReturn(0);
}

/*@
  DMPlexCreateProcessSF - Create an SF which just has process connectivity

  Collective on dm

  Input Parameters:
+ dm      - The DM
- sfPoint - The PetscSF which encodes point connectivity

  Output Parameters:
+ processRanks - A list of process neighbors, or NULL
- sfProcess    - An SF encoding the process connectivity, or NULL

  Level: developer

.seealso: PetscSFCreate(), DMPlexCreateTwoSidedProcessSF()
@*/
PetscErrorCode DMPlexCreateProcessSF(DM dm, PetscSF sfPoint, IS *processRanks, PetscSF *sfProcess)
{
  PetscInt           numRoots, numLeaves, l;
  const PetscInt    *localPoints;
  const PetscSFNode *remotePoints;
  PetscInt          *localPointsNew;
  PetscSFNode       *remotePointsNew;
  PetscInt          *ranks, *ranksNew;
  PetscMPIInt        size;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidHeaderSpecific(sfPoint, PETSCSF_CLASSID, 2);
  if (processRanks) {PetscValidPointer(processRanks, 3);}
  if (sfProcess)    {PetscValidPointer(sfProcess, 4);}
  ierr = MPI_Comm_size(PetscObjectComm((PetscObject) dm), &size);CHKERRQ(ierr);
  ierr = PetscSFGetGraph(sfPoint, &numRoots, &numLeaves, &localPoints, &remotePoints);CHKERRQ(ierr);
  ierr = PetscMalloc1(numLeaves, &ranks);CHKERRQ(ierr);
  for (l = 0; l < numLeaves; ++l) {
    ranks[l] = remotePoints[l].rank;
  }
  ierr = PetscSortRemoveDupsInt(&numLeaves, ranks);CHKERRQ(ierr);
  ierr = PetscMalloc1(numLeaves, &ranksNew);CHKERRQ(ierr);
  ierr = PetscMalloc1(numLeaves, &localPointsNew);CHKERRQ(ierr);
  ierr = PetscMalloc1(numLeaves, &remotePointsNew);CHKERRQ(ierr);
  for (l = 0; l < numLeaves; ++l) {
    ranksNew[l]              = ranks[l];
    localPointsNew[l]        = l;
    remotePointsNew[l].index = 0;
    remotePointsNew[l].rank  = ranksNew[l];
  }
  ierr = PetscFree(ranks);CHKERRQ(ierr);
  if (processRanks) {ierr = ISCreateGeneral(PetscObjectComm((PetscObject)dm), numLeaves, ranksNew, PETSC_OWN_POINTER, processRanks);CHKERRQ(ierr);}
  else              {ierr = PetscFree(ranksNew);CHKERRQ(ierr);}
  if (sfProcess) {
    ierr = PetscSFCreate(PetscObjectComm((PetscObject)dm), sfProcess);CHKERRQ(ierr);
    ierr = PetscObjectSetName((PetscObject) *sfProcess, "Process SF");CHKERRQ(ierr);
    ierr = PetscSFSetFromOptions(*sfProcess);CHKERRQ(ierr);
    ierr = PetscSFSetGraph(*sfProcess, size, numLeaves, localPointsNew, PETSC_OWN_POINTER, remotePointsNew, PETSC_OWN_POINTER);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerCreateSF(DMPlexCellRefiner cr, DM rdm)
{
  DM                 dm = cr->dm;
  DMPlexCellRefiner *crRem;
  PetscSF            sf, sfNew, sfProcess;
  IS                 processRanks;
  MPI_Datatype       ctType;
  PetscInt           numRoots, numLeaves, numLeavesNew = 0, l, m;
  const PetscInt    *localPoints, *neighbors;
  const PetscSFNode *remotePoints;
  PetscInt          *localPointsNew;
  PetscSFNode       *remotePointsNew;
  PetscInt          *ctStartRem, *ctStartNewRem;
  PetscInt           ctSize = DM_NUM_POLYTOPES+1, numNeighbors, n, pStartNew, pEndNew, pNew, pNewRem;
  PetscErrorCode     ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetChart(rdm, &pStartNew, &pEndNew);CHKERRQ(ierr);
  ierr = DMGetPointSF(dm, &sf);CHKERRQ(ierr);
  ierr = DMGetPointSF(rdm, &sfNew);CHKERRQ(ierr);
  /* Calculate size of new SF */
  ierr = PetscSFGetGraph(sf, &numRoots, &numLeaves, &localPoints, &remotePoints);CHKERRQ(ierr);
  if (numRoots < 0) PetscFunctionReturn(0);
  for (l = 0; l < numLeaves; ++l) {
    const PetscInt  p = localPoints[l];
    DMPolytopeType  ct;
    DMPolytopeType *rct;
    PetscInt       *rsize, *rcone, *rornt;
    PetscInt        Nct, n;

    ierr = DMPlexGetCellType(dm, p, &ct);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
    for (n = 0; n < Nct; ++n) numLeavesNew += rsize[n];
  }
  /* Communicate ctStart and cStartNew for each remote rank */
  ierr = DMPlexCreateProcessSF(dm, sf, &processRanks, &sfProcess);CHKERRQ(ierr);
  ierr = ISGetLocalSize(processRanks, &numNeighbors);CHKERRQ(ierr);
  ierr = PetscMalloc2(ctSize*numNeighbors, &ctStartRem, ctSize*numNeighbors, &ctStartNewRem);CHKERRQ(ierr);
  ierr = MPI_Type_contiguous(ctSize, MPIU_INT, &ctType);CHKERRQ(ierr);
  ierr = MPI_Type_commit(&ctType);CHKERRQ(ierr);
  ierr = PetscSFBcastBegin(sfProcess, ctType, cr->ctStart, ctStartRem);CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfProcess, ctType, cr->ctStart, ctStartRem);CHKERRQ(ierr);
  ierr = PetscSFBcastBegin(sfProcess, ctType, cr->ctStartNew, ctStartNewRem);CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfProcess, ctType, cr->ctStartNew, ctStartNewRem);CHKERRQ(ierr);
  ierr = MPI_Type_free(&ctType);CHKERRQ(ierr);
  ierr = PetscSFDestroy(&sfProcess);CHKERRQ(ierr);
  ierr = PetscMalloc1(numNeighbors, &crRem);CHKERRQ(ierr);
  for (n = 0; n < numNeighbors; ++n) {
    ierr = DMPlexCellRefinerCreate(dm, &crRem[n]);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerSetStarts(crRem[n], &ctStartRem[n*ctSize], &ctStartNewRem[n*ctSize]);
    ierr = DMPlexCellRefinerSetUp(crRem[n]);CHKERRQ(ierr);
  }
  ierr = PetscFree2(ctStartRem, ctStartNewRem);CHKERRQ(ierr);
  /* Calculate new point SF */
  ierr = PetscMalloc1(numLeavesNew, &localPointsNew);CHKERRQ(ierr);
  ierr = PetscMalloc1(numLeavesNew, &remotePointsNew);CHKERRQ(ierr);
  ierr = ISGetIndices(processRanks, &neighbors);CHKERRQ(ierr);
  for (l = 0, m = 0; l < numLeaves; ++l) {
    PetscInt        p       = localPoints[l];
    PetscInt        pRem    = remotePoints[l].index;
    PetscMPIInt     rankRem = remotePoints[l].rank;
    DMPolytopeType  ct;
    DMPolytopeType *rct;
    PetscInt       *rsize, *rcone, *rornt;
    PetscInt        neighbor, Nct, n, r;

    ierr = PetscFindInt(rankRem, numNeighbors, neighbors, &neighbor);CHKERRQ(ierr);
    if (neighbor < 0) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Could not locate remote rank %D", rankRem);
    ierr = DMPlexGetCellType(dm, p, &ct);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
    for (n = 0; n < Nct; ++n) {
      for (r = 0; r < rsize[n]; ++r) {
        ierr = DMPlexCellRefinerGetNewPoint(cr, ct, rct[n], p, r, &pNew);CHKERRQ(ierr);
        ierr = DMPlexCellRefinerGetNewPoint(crRem[neighbor], ct, rct[n], pRem, r, &pNewRem);CHKERRQ(ierr);
        localPointsNew[m]        = pNew;
        remotePointsNew[m].index = pNewRem;
        remotePointsNew[m].rank  = rankRem;
        ++m;
      }
    }
  }
  for (n = 0; n < numNeighbors; ++n) {ierr = DMPlexCellRefinerDestroy(&crRem[n]);CHKERRQ(ierr);}
  ierr = PetscFree(crRem);CHKERRQ(ierr);
  if (m != numLeavesNew) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Number of leaf point %D should be %D", m, numLeavesNew);
  ierr = ISRestoreIndices(processRanks, &neighbors);CHKERRQ(ierr);
  ierr = ISDestroy(&processRanks);CHKERRQ(ierr);
  {
    PetscSFNode *rp, *rtmp;
    PetscInt    *lp, *idx, *ltmp, i;

    /* SF needs sorted leaves to correct calculate Gather */
    ierr = PetscMalloc1(numLeavesNew, &idx);CHKERRQ(ierr);
    ierr = PetscMalloc1(numLeavesNew, &lp);CHKERRQ(ierr);
    ierr = PetscMalloc1(numLeavesNew, &rp);CHKERRQ(ierr);
    for (i = 0; i < numLeavesNew; ++i) {
      if ((localPointsNew[i] < pStartNew) || (localPointsNew[i] >= pEndNew)) SETERRQ4(PETSC_COMM_SELF, PETSC_ERR_ARG_OUTOFRANGE, "Local SF point %D (%D) not in [%D, %D)", localPointsNew[i], i, pStartNew, pEndNew);
      idx[i] = i;
    }
    ierr = PetscSortIntWithPermutation(numLeavesNew, localPointsNew, idx);CHKERRQ(ierr);
    for (i = 0; i < numLeavesNew; ++i) {
      lp[i] = localPointsNew[idx[i]];
      rp[i] = remotePointsNew[idx[i]];
    }
    ltmp            = localPointsNew;
    localPointsNew  = lp;
    rtmp            = remotePointsNew;
    remotePointsNew = rp;
    ierr = PetscFree(idx);CHKERRQ(ierr);
    ierr = PetscFree(ltmp);CHKERRQ(ierr);
    ierr = PetscFree(rtmp);CHKERRQ(ierr);
  }
  ierr = PetscSFSetGraph(sfNew, pEndNew-pStartNew, numLeavesNew, localPointsNew, PETSC_OWN_POINTER, remotePointsNew, PETSC_OWN_POINTER);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMPlexCellRefinerCreateLabels(DMPlexCellRefiner cr, DM rdm)
{
  DM             dm = cr->dm;
  PetscInt       numLabels, l;
  PetscInt       pNew;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetNumLabels(dm, &numLabels);CHKERRQ(ierr);
  for (l = 0; l < numLabels; ++l) {
    DMLabel         label, labelNew;
    const char     *lname;
    PetscBool       isDepth, isCellType;
    IS              valueIS;
    const PetscInt *values;
    PetscInt        defVal;
    PetscInt        numValues, val;

    ierr = DMGetLabelName(dm, l, &lname);CHKERRQ(ierr);
    ierr = PetscStrcmp(lname, "depth", &isDepth);CHKERRQ(ierr);
    if (isDepth) continue;
    ierr = PetscStrcmp(lname, "celltype", &isCellType);CHKERRQ(ierr);
    if (isCellType) continue;
    ierr = DMCreateLabel(rdm, lname);CHKERRQ(ierr);
    ierr = DMGetLabel(dm, lname, &label);CHKERRQ(ierr);
    ierr = DMGetLabel(rdm, lname, &labelNew);CHKERRQ(ierr);
    ierr = DMLabelGetDefaultValue(label, &defVal);CHKERRQ(ierr);
    ierr = DMLabelSetDefaultValue(labelNew, defVal);CHKERRQ(ierr);
    ierr = DMLabelGetValueIS(label, &valueIS);CHKERRQ(ierr);
    ierr = ISGetLocalSize(valueIS, &numValues);CHKERRQ(ierr);
    ierr = ISGetIndices(valueIS, &values);CHKERRQ(ierr);
    for (val = 0; val < numValues; ++val) {
      IS              pointIS;
      const PetscInt *points;
      PetscInt        numPoints, p;

      /* Ensure refined label is created with same number of strata as
       * original (even if no entries here). */
      ierr = DMLabelAddStratum(labelNew, values[val]);CHKERRQ(ierr);
      ierr = DMLabelGetStratumIS(label, values[val], &pointIS);CHKERRQ(ierr);
      ierr = ISGetLocalSize(pointIS, &numPoints);CHKERRQ(ierr);
      ierr = ISGetIndices(pointIS, &points);CHKERRQ(ierr);
      for (p = 0; p < numPoints; ++p) {
        const PetscInt  point = points[p];
        DMPolytopeType  ct;
        DMPolytopeType *rct;
        PetscInt       *rsize, *rcone, *rornt;
        PetscInt        Nct, n, r;

        ierr = DMPlexGetCellType(dm, point, &ct);CHKERRQ(ierr);
        ierr = DMPlexCellRefinerRefine(cr, ct, &Nct, &rct, &rsize, &rcone, &rornt);CHKERRQ(ierr);
        for (n = 0; n < Nct; ++n) {
          for (r = 0; r < rsize[n]; ++r) {
            ierr = DMPlexCellRefinerGetNewPoint(cr, ct, rct[n], point, r, &pNew);CHKERRQ(ierr);
            ierr = DMLabelSetValue(labelNew, pNew, values[val]);CHKERRQ(ierr);
          }
        }
      }
      ierr = ISRestoreIndices(pointIS, &points);CHKERRQ(ierr);
      ierr = ISDestroy(&pointIS);CHKERRQ(ierr);
    }
    ierr = ISRestoreIndices(valueIS, &values);CHKERRQ(ierr);
    ierr = ISDestroy(&valueIS);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* This will only work for interpolated meshes */
PetscErrorCode DMPlexRefineUniform(DM dm, DMPlexCellRefiner cr, DM *dmRefined)
{
  DM              rdm;
  PetscInt        dim, embedDim, depth;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  PetscValidHeader(cr, 1);
  ierr = DMCreate(PetscObjectComm((PetscObject)dm), &rdm);CHKERRQ(ierr);
  ierr = DMSetType(rdm, DMPLEX);CHKERRQ(ierr);
  ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMSetDimension(rdm, dim);CHKERRQ(ierr);
  ierr = DMGetCoordinateDim(dm, &embedDim);CHKERRQ(ierr);
  ierr = DMSetCoordinateDim(rdm, embedDim);CHKERRQ(ierr);
  /* Calculate number of new points of each depth */
  ierr = DMPlexGetDepth(dm, &depth);CHKERRQ(ierr);
  if (depth >= 0 && dim != depth) SETERRQ(PetscObjectComm((PetscObject) dm), PETSC_ERR_ARG_WRONG, "Mesh must be interpolated for regular refinement");
  /* Step 1: Set chart */
  ierr = DMPlexSetChart(rdm, 0, cr->ctStartNew[cr->ctOrder[DM_NUM_POLYTOPES]]);CHKERRQ(ierr);
  /* Step 2: Set cone/support sizes (automatically stratifies) */
  ierr = DMPlexCellRefinerSetConeSizes(cr, rdm);CHKERRQ(ierr);
  /* Step 3: Setup refined DM */
  ierr = DMSetUp(rdm);CHKERRQ(ierr);
  /* Step 4: Set cones and supports (automatically symmetrizes) */
  ierr = DMPlexCellRefinerSetCones(cr, rdm);CHKERRQ(ierr);
  /* Step 5: Create pointSF */
  ierr = DMPlexCellRefinerCreateSF(cr, rdm);CHKERRQ(ierr);
  /* Step 6: Create labels */
  ierr = DMPlexCellRefinerCreateLabels(cr, rdm);CHKERRQ(ierr);
  /* Step 7: Set coordinates */
  ierr = DMPlexCellRefinerSetCoordinates(cr, rdm);CHKERRQ(ierr);
  *dmRefined = rdm;
  PetscFunctionReturn(0);
}

/*@
  DMPlexCreateCoarsePointIS - Creates an IS covering the coarse DM chart with the fine points as data

  Input Parameter:
. dm - The coarse DM

  Output Parameter:
. fpointIS - The IS of all the fine points which exist in the original coarse mesh

  Level: developer

.seealso: DMRefine(), DMPlexSetRefinementUniform(), DMPlexCreateSubpointIS()
@*/
PetscErrorCode DMPlexCreateCoarsePointIS(DM dm, IS *fpointIS)
{
  DMPlexCellRefiner cr;
  PetscInt         *fpoints;
  PetscInt          pStart, pEnd, p, vStart, vEnd, v;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetChart(dm, &pStart, &pEnd);CHKERRQ(ierr);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  ierr = DMPlexCellRefinerCreate(dm, &cr);CHKERRQ(ierr);
  ierr = DMPlexCellRefinerSetUp(cr);CHKERRQ(ierr);
  ierr = PetscMalloc1(pEnd-pStart, &fpoints);CHKERRQ(ierr);
  for (p = 0; p < pEnd-pStart; ++p) fpoints[p] = -1;
  for (v = vStart; v < vEnd; ++v) {
    PetscInt vNew = -1; /* silent overzelous may be used uninitialized */

    ierr = DMPlexCellRefinerGetNewPoint(cr, DM_POLYTOPE_POINT, DM_POLYTOPE_POINT, p, 0, &vNew);CHKERRQ(ierr);
    fpoints[v-pStart] = vNew;
  }
  ierr = DMPlexCellRefinerDestroy(&cr);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF, pEnd-pStart, fpoints, PETSC_OWN_POINTER, fpointIS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
  DMPlexSetRefinementUniform - Set the flag for uniform refinement

  Input Parameters:
+ dm - The DM
- refinementUniform - The flag for uniform refinement

  Level: developer

.seealso: DMRefine(), DMPlexGetRefinementUniform(), DMPlexGetRefinementLimit(), DMPlexSetRefinementLimit()
@*/
PetscErrorCode DMPlexSetRefinementUniform(DM dm, PetscBool refinementUniform)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  mesh->refinementUniform = refinementUniform;
  PetscFunctionReturn(0);
}

/*@
  DMPlexGetRefinementUniform - Retrieve the flag for uniform refinement

  Input Parameter:
. dm - The DM

  Output Parameter:
. refinementUniform - The flag for uniform refinement

  Level: developer

.seealso: DMRefine(), DMPlexSetRefinementUniform(), DMPlexGetRefinementLimit(), DMPlexSetRefinementLimit()
@*/
PetscErrorCode DMPlexGetRefinementUniform(DM dm, PetscBool *refinementUniform)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidPointer(refinementUniform,  2);
  *refinementUniform = mesh->refinementUniform;
  PetscFunctionReturn(0);
}

/*@
  DMPlexSetRefinementLimit - Set the maximum cell volume for refinement

  Input Parameters:
+ dm - The DM
- refinementLimit - The maximum cell volume in the refined mesh

  Level: developer

.seealso: DMRefine(), DMPlexGetRefinementLimit(), DMPlexGetRefinementUniform(), DMPlexSetRefinementUniform()
@*/
PetscErrorCode DMPlexSetRefinementLimit(DM dm, PetscReal refinementLimit)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  mesh->refinementLimit = refinementLimit;
  PetscFunctionReturn(0);
}

/*@
  DMPlexGetRefinementLimit - Retrieve the maximum cell volume for refinement

  Input Parameter:
. dm - The DM

  Output Parameter:
. refinementLimit - The maximum cell volume in the refined mesh

  Level: developer

.seealso: DMRefine(), DMPlexSetRefinementLimit(), DMPlexGetRefinementUniform(), DMPlexSetRefinementUniform()
@*/
PetscErrorCode DMPlexGetRefinementLimit(DM dm, PetscReal *refinementLimit)
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidPointer(refinementLimit,  2);
  /* if (mesh->refinementLimit < 0) = getMaxVolume()/2.0; */
  *refinementLimit = mesh->refinementLimit;
  PetscFunctionReturn(0);
}

/*@
  DMPlexSetRefinementFunction - Set the function giving the maximum cell volume for refinement

  Input Parameters:
+ dm - The DM
- refinementFunc - Function giving the maximum cell volume in the refined mesh

  Note: The calling sequence is refinementFunc(coords, limit)
$ coords - Coordinates of the current point, usually a cell centroid
$ limit  - The maximum cell volume for a cell containing this point

  Level: developer

.seealso: DMRefine(), DMPlexGetRefinementFunction(), DMPlexGetRefinementUniform(), DMPlexSetRefinementUniform(), DMPlexGetRefinementLimit(), DMPlexSetRefinementLimit()
@*/
PetscErrorCode DMPlexSetRefinementFunction(DM dm, PetscErrorCode (*refinementFunc)(const PetscReal [], PetscReal *))
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  mesh->refinementFunc = refinementFunc;
  PetscFunctionReturn(0);
}

/*@
  DMPlexGetRefinementFunction - Get the function giving the maximum cell volume for refinement

  Input Parameter:
. dm - The DM

  Output Parameter:
. refinementFunc - Function giving the maximum cell volume in the refined mesh

  Note: The calling sequence is refinementFunc(coords, limit)
$ coords - Coordinates of the current point, usually a cell centroid
$ limit  - The maximum cell volume for a cell containing this point

  Level: developer

.seealso: DMRefine(), DMPlexSetRefinementFunction(), DMPlexGetRefinementUniform(), DMPlexSetRefinementUniform(), DMPlexGetRefinementLimit(), DMPlexSetRefinementLimit()
@*/
PetscErrorCode DMPlexGetRefinementFunction(DM dm, PetscErrorCode (**refinementFunc)(const PetscReal [], PetscReal *))
{
  DM_Plex *mesh = (DM_Plex*) dm->data;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(dm, DM_CLASSID, 1);
  PetscValidPointer(refinementFunc,  2);
  *refinementFunc = mesh->refinementFunc;
  PetscFunctionReturn(0);
}

PetscErrorCode DMRefine_Plex(DM dm, MPI_Comm comm, DM *dmRefined)
{
  PetscBool         isUniform;
  DMPlexCellRefiner cr;
  PetscErrorCode    ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetRefinementUniform(dm, &isUniform);CHKERRQ(ierr);
  ierr = DMViewFromOptions(dm, NULL, "-initref_dm_view");CHKERRQ(ierr);
  if (isUniform) {
    PetscBool localized;

    ierr = DMPlexCellRefinerCreate(dm, &cr);CHKERRQ(ierr);
    ierr = DMPlexCellRefinerSetUp(cr);CHKERRQ(ierr);
    ierr = DMGetCoordinatesLocalized(dm, &localized);CHKERRQ(ierr);
    ierr = DMPlexRefineUniform(dm, cr, dmRefined);CHKERRQ(ierr);
    ierr = DMPlexSetRegularRefinement(*dmRefined, PETSC_TRUE);CHKERRQ(ierr);
    ierr = DMCopyBoundary(dm, *dmRefined);CHKERRQ(ierr);
    if (localized) {ierr = DMLocalizeCoordinates(*dmRefined);CHKERRQ(ierr);}
    ierr = DMPlexCellRefinerDestroy(&cr);CHKERRQ(ierr);
  } else {
    ierr = DMPlexRefine_Internal(dm, NULL, dmRefined);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode DMRefineHierarchy_Plex(DM dm, PetscInt nlevels, DM dmRefined[])
{
  DM             cdm = dm;
  PetscInt       r;
  PetscBool      isUniform, localized;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMPlexGetRefinementUniform(dm, &isUniform);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocalized(dm, &localized);CHKERRQ(ierr);
  if (isUniform) {
    for (r = 0; r < nlevels; ++r) {
      DMPlexCellRefiner cr;

      ierr = DMPlexCellRefinerCreate(cdm, &cr);CHKERRQ(ierr);
      ierr = DMPlexCellRefinerSetUp(cr);CHKERRQ(ierr);
      ierr = DMPlexRefineUniform(cdm, cr, &dmRefined[r]);CHKERRQ(ierr);
      ierr = DMSetCoarsenLevel(dmRefined[r], cdm->leveldown);CHKERRQ(ierr);
      ierr = DMSetRefineLevel(dmRefined[r], cdm->levelup+1);CHKERRQ(ierr);
      ierr = DMCopyBoundary(cdm, dmRefined[r]);CHKERRQ(ierr);
      if (localized) {ierr = DMLocalizeCoordinates(dmRefined[r]);CHKERRQ(ierr);}
      ierr = DMSetCoarseDM(dmRefined[r], cdm);CHKERRQ(ierr);
      ierr = DMPlexSetRegularRefinement(dmRefined[r], PETSC_TRUE);CHKERRQ(ierr);
      cdm  = dmRefined[r];
      ierr = DMPlexCellRefinerDestroy(&cr);CHKERRQ(ierr);
    }
  } else {
    for (r = 0; r < nlevels; ++r) {
      ierr = DMRefine(cdm, PetscObjectComm((PetscObject) dm), &dmRefined[r]);CHKERRQ(ierr);
      ierr = DMCopyBoundary(cdm, dmRefined[r]);CHKERRQ(ierr);
      if (localized) {ierr = DMLocalizeCoordinates(dmRefined[r]);CHKERRQ(ierr);}
      ierr = DMSetCoarseDM(dmRefined[r], cdm);CHKERRQ(ierr);
      cdm  = dmRefined[r];
    }
  }
  PetscFunctionReturn(0);
}
