#include <petsc/private/fortranimpl.h>


#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define petscinfo_ PETSCINFO
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE) && !defined(FORTRANDOUBLEUNDERSCORE)
#define petscinfo_ petscinfo
#endif

static PetscErrorCode PetscFixSlashN(const char *in, char **out)
{
  PetscErrorCode ierr;
  PetscInt       i;
  size_t         len;

  PetscFunctionBegin;
  ierr = PetscStrallocpy(in,out);CHKERRQ(ierr);
  ierr = PetscStrlen(*out,&len);CHKERRQ(ierr);
  for (i=0; i<(int)len-1; i++) {
    if ((*out)[i] == '\\' && (*out)[i+1] == 'n') {(*out)[i] = ' '; (*out)[i+1] = '\n';}
  }
  PetscFunctionReturn(0);
}

PETSC_EXTERN void petscinfosetfile_(char* filename, char* mode,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len1,PETSC_FORTRAN_CHARLEN_T len2)
{
  char *t1, *t2;

  FIXCHAR(filename,len1,t1);
  FIXCHAR(mode,len2,t2);
  *ierr = PetscInfoSetFile(t1,t2); if (*ierr) return;
  FREECHAR(filename,t1);
  FREECHAR(mode,t2);
}

PETSC_EXTERN void petscinfogetclass_(char *classname,PetscBool **found,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(classname,len,t);
  *ierr = PetscInfoGetClass(t,*found);if (*ierr) return;
  FREECHAR(classname,t);
}

PETSC_EXTERN void petscinfoprocessclass_(char *classname,PetscInt *numClassID,PetscClassId *classIDs[],PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len)
{
  char *t;

  FIXCHAR(classname,len,t);
  *ierr = PetscInfoProcessClass(t,*numClassID,*classIDs); if (*ierr) return;
  FREECHAR(classname,t);
}

PETSC_EXTERN void petscinfo_(char* text,PetscErrorCode *ierr,PETSC_FORTRAN_CHARLEN_T len1)
{
  char *c1,*tmp;

  FIXCHAR(text,len1,c1);
  *ierr = PetscFixSlashN(c1,&tmp);if (*ierr) return;
  FREECHAR(text,c1);
  *ierr = PetscInfo(NULL,tmp);if (*ierr) return;
  *ierr = PetscFree(tmp);
}
