#if !defined(PETSCSECTIONTYPES_H)
#define PETSCSECTIONTYPES_H

/*S
  PetscSection - Mapping from integers in a designated range to contiguous sets of integers.

  In contrast to IS, which maps from integers to single integers, the range of a PetscSection is in the space of
  contiguous sets of integers. These ranges are frequently interpreted as domains of other array-like objects,
  especially other PetscSections, Vecs, and ISs. The domain is set with PetscSectionSetChart() and does not need to
  start at 0. For each point in the domain of a PetscSection, the output set is represented through an offset and a
  count, which are set using PetscSectionSetOffset() and PetscSectionSetDof() respectively. Lookup is typically using
  accessors or routines like VecGetValuesSection().

  Level: beginner

.seealso:  PetscSectionCreate(), PetscSectionDestroy(), PetscSectionSym
S*/
typedef struct _p_PetscSection *PetscSection;

/*S
  PetscSectionSym - Symmetries of the data referenced by a PetscSection.

  Often the order of data index by a PetscSection is meaningful, and describes additional structure, such as points on a
  line, grid, or lattice.  If the data is accessed from a different "orientation", then the image of the data under
  access then undergoes a symmetry transformation.  A PetscSectionSym specifies these symmetries.  The types of
  symmetries that can be specified are of the form R * P, where R is a diagonal matrix of scalars, and P is a permutation.

  Level: developer

.seealso: PetscSection, PetscSectionSymCreate(), PetscSectionSymDestroy(), PetscSectionSetSym(), PetscSectionGetSym(), PetscSectionSetFieldSym(), 
          PetscSectionGetFieldSym(), PetscSectionGetSymPoints(), PetscSectionSymType, PetscSectionSymSetType(), PetscSectionSymGetType()
S*/
typedef struct _p_PetscSectionSym *PetscSectionSym;

/*J
  PetscSectionSymType - String with the name of a PetscSectionSym type.

  Level: developer

  Notes:
    PetscSectionSym has no default implementation, but is used by DM in PetscSectionSymCreateLabel().

.seealso: PetscSectionSymSetType(), PetscSectionSymGetType(), PetscSectionSym, PetscSectionSymCreate(), PetscSectionSymRegister()
J*/
typedef const char *PetscSectionSymType;

#endif
