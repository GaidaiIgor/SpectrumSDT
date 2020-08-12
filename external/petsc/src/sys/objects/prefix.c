
/*
     Provides utility routines for manulating any type of PETSc object.
*/
#include <petsc/private/petscimpl.h>  /*I   "petscsys.h"    I*/

/*@C
   PetscObjectGetOptions - Gets the options database used by the object. Call immediately after creating the object.

   Collective on PetscObject

   Input Parameter:
.  obj - any PETSc object, for example a Vec, Mat or KSP.

   Output Parameter:
.  options - the options database

   Notes:
    if this is not called the object will use the default options database

  Level: advanced

.seealso: PetscOptionsCreate(), PetscOptionsDestroy(), PetscObjectSetOptionsPrefix(), PetscObjectAppendOptionsPrefix(), PetscObjectPrependOptionsPrefix(),
          PetscObjectGetOptionsPrefix(), PetscObjectSetOptions()

@*/
PetscErrorCode  PetscObjectGetOptions(PetscObject obj,PetscOptions *options)
{
  PetscFunctionBegin;
  PetscValidHeader(obj,1);
  *options = obj->options;
  PetscFunctionReturn(0);
}

/*@C
   PetscObjectSetOptions - Sets the options database used by the object. Call immediately after creating the object.

   Collective on PetscObject

   Input Parameters:
+  obj - any PETSc object, for example a Vec, Mat or KSP.
-  options - the options database, use NULL for default

   Notes:
    if this is not called the object will use the default options database

  Level: advanced

.seealso: PetscOptionsCreate(), PetscOptionsDestroy(), PetscObjectSetOptionsPrefix(), PetscObjectAppendOptionsPrefix(), PetscObjectPrependOptionsPrefix(),
          PetscObjectGetOptionsPrefix(), PetscObjectGetOptions()

@*/
PetscErrorCode  PetscObjectSetOptions(PetscObject obj,PetscOptions options)
{
  PetscFunctionBegin;
  PetscValidHeader(obj,1);
  obj->options = options;
  PetscFunctionReturn(0);
}

/*@C
   PetscObjectSetOptionsPrefix - Sets the prefix used for searching for all
   options of PetscObjectType in the database.

   Collective on Object

   Input Parameters:
+  obj - any PETSc object, for example a Vec, Mat or KSP.
-  prefix - the prefix string to prepend to option requests of the object.

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

  Level: advanced

.seealso: PetscOptionsCreate(), PetscOptionsDestroy(), PetscObjectAppendOptionsPrefix(), PetscObjectPrependOptionsPrefix(),
          PetscObjectGetOptionsPrefix(), TSSetOptionsPrefix(), SNESSetOptionsPrefix(), KSPSetOptionsPrefix()

@*/
PetscErrorCode  PetscObjectSetOptionsPrefix(PetscObject obj,const char prefix[])
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeader(obj,1);
  if (!prefix) {
    ierr = PetscFree(obj->prefix);CHKERRQ(ierr);
  } else {
    if (prefix[0] == '-') SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Options prefix should not begin with a hypen");
    if (prefix != obj->prefix) {
      ierr = PetscFree(obj->prefix);CHKERRQ(ierr);
      ierr = PetscStrallocpy(prefix,&obj->prefix);CHKERRQ(ierr);
    }
  }
  PetscFunctionReturn(0);
}

/*@C
   PetscObjectAppendOptionsPrefix - Sets the prefix used for searching for all
   options of PetscObjectType in the database.

   Input Parameters:
+  obj - any PETSc object, for example a Vec, Mat or KSP.
-  prefix - the prefix string to prepend to option requests of the object.

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

  Level: advanced

.seealso: PetscOptionsCreate(), PetscOptionsDestroy(), PetscObjectSetOptionsPrefix(), PetscObjectPrependOptionsPrefix(),
          PetscObjectGetOptionsPrefix(), TSAppendOptionsPrefix(), SNESAppendOptionsPrefix(), KSPAppendOptionsPrefix()

@*/
PetscErrorCode  PetscObjectAppendOptionsPrefix(PetscObject obj,const char prefix[])
{
  char           *buf = obj->prefix;
  PetscErrorCode ierr;
  size_t         len1,len2;

  PetscFunctionBegin;
  PetscValidHeader(obj,1);
  if (!prefix) PetscFunctionReturn(0);
  if (!buf) {
    ierr = PetscObjectSetOptionsPrefix(obj,prefix);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (prefix[0] == '-') SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Options prefix should not begin with a hypen");

  ierr = PetscStrlen(prefix,&len1);CHKERRQ(ierr);
  ierr = PetscStrlen(buf,&len2);CHKERRQ(ierr);
  ierr = PetscMalloc1(1+len1+len2,&obj->prefix);CHKERRQ(ierr);
  ierr = PetscStrcpy(obj->prefix,buf);CHKERRQ(ierr);
  ierr = PetscStrcat(obj->prefix,prefix);CHKERRQ(ierr);
  ierr = PetscFree(buf);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PetscObjectGetOptionsPrefix - Gets the prefix of the PetscObject.

   Input Parameters:
.  obj - any PETSc object, for example a Vec, Mat or KSP.

   Output Parameters:
.  prefix - pointer to the prefix string used is returned

  Level: advanced

.seealso: PetscOptionsCreate(), PetscOptionsDestroy(), PetscObjectSetOptionsPrefix(), PetscObjectAppendOptionsPrefix(), PetscObjectPrependOptionsPrefix(),
          TSGetOptionsPrefix(), SNESGetOptionsPrefix(), KSPGetOptionsPrefix()

@*/
PetscErrorCode  PetscObjectGetOptionsPrefix(PetscObject obj,const char *prefix[])
{
  PetscFunctionBegin;
  PetscValidHeader(obj,1);
  PetscValidPointer(prefix,2);
  *prefix = obj->prefix;
  PetscFunctionReturn(0);
}

/*@C
   PetscObjectPrependOptionsPrefix - Sets the prefix used for searching for all
   options of PetscObjectType in the database.

   Input Parameters:
+  obj - any PETSc object, for example a Vec, Mat or KSP.
-  prefix - the prefix string to prepend to option requests of the object.

   Notes:
   A hyphen (-) must NOT be given at the beginning of the prefix name.
   The first character of all runtime options is AUTOMATICALLY the
   hyphen.

  Level: advanced

.seealso: PetscOptionsCreate(), PetscOptionsDestroy(), PetscObjectSetOptionsPrefix(), PetscObjectAppendOptionsPrefix(),
          PetscObjectGetOptionsPrefix(), TSPrependOptionsPrefix(), SNESPrependOptionsPrefix(), KSPPrependOptionsPrefix()

@*/
PetscErrorCode  PetscObjectPrependOptionsPrefix(PetscObject obj,const char prefix[])
{
  char           *buf;
  size_t         len1,len2;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeader(obj,1);
  buf = obj->prefix;
  if (!prefix) PetscFunctionReturn(0);
  if (!buf) {
    ierr = PetscObjectSetOptionsPrefix(obj,prefix);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  if (prefix[0] == '-') SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"Options prefix should not begin with a hypen");

  ierr = PetscStrlen(prefix,&len1);CHKERRQ(ierr);
  ierr = PetscStrlen(buf,&len2);CHKERRQ(ierr);
  ierr = PetscMalloc1(1+len1+len2,&obj->prefix);CHKERRQ(ierr);
  ierr = PetscStrcpy(obj->prefix,prefix);CHKERRQ(ierr);
  ierr = PetscStrcat(obj->prefix,buf);CHKERRQ(ierr);
  ierr = PetscFree(buf);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

