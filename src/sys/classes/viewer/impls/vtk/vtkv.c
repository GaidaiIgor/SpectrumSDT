#include <../src/sys/classes/viewer/impls/vtk/vtkvimpl.h> /*I "petscviewer.h" I*/

/*MC
    PetscViewerVTKWriteFunction - functional form used to provide writer to the PetscViewerVTK

     Synopsis:
     #include <petscviewer.h>
     PetscViewerVTKWriteFunction(PetscObject object,PetscViewer viewer)

     Input Parameters:
+      object - the PETSc object to be written
-      viewer - viewer it is to be written to

   Level: developer

.seealso:   PetscViewerVTKAddField()
M*/

/*@C
   PetscViewerVTKAddField - Add a field to the viewer

   Collective

   Input Arguments:
+ viewer - VTK viewer
. dm - DM on which Vec lives
. PetscViewerVTKWriteFunction - function to write this Vec
. fieldnum - which field of the DM to write (PETSC_DEFAULT if the whle vector should be written)
. fieldtype - Either PETSC_VTK_POINT_FIELD or PETSC_VTK_CELL_FIELD
. checkdm - whether to check for identical dm arguments as fields are added
- vec - Vec from which to write

   Note:
   This routine keeps exclusive ownership of the Vec. The caller should not use or destroy the Vec after adding it.

   Level: developer

.seealso: PetscViewerVTKOpen(), DMDAVTKWriteAll(), PetscViewerVTKWriteFunction, PetscViewerVTKGetDM()
@*/
PetscErrorCode PetscViewerVTKAddField(PetscViewer viewer,PetscObject dm,PetscErrorCode (*PetscViewerVTKWriteFunction)(PetscObject,PetscViewer),PetscInt fieldnum,PetscViewerVTKFieldType fieldtype,PetscBool checkdm,PetscObject vec)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,1);
  PetscValidHeader(dm,2);
  PetscValidHeader(vec,4);
  ierr = PetscUseMethod(viewer,"PetscViewerVTKAddField_C",(PetscViewer,PetscObject,PetscErrorCode (*)(PetscObject,PetscViewer),PetscInt,PetscViewerVTKFieldType,PetscBool,PetscObject),(viewer,dm,PetscViewerVTKWriteFunction,fieldnum,fieldtype,checkdm,vec));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PetscViewerVTKGetDM - get the DM associated with the viewer

   Collective

   Input Arguments:
+ viewer - VTK viewer
- dm - DM associated with the viewer (as PetscObject)

   Level: developer

.seealso: PetscViewerVTKOpen(), DMDAVTKWriteAll(), PetscViewerVTKWriteFunction, PetscViewerVTKAddField()
@*/
PetscErrorCode PetscViewerVTKGetDM(PetscViewer viewer,PetscObject *dm)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidHeaderSpecific(viewer,PETSC_VIEWER_CLASSID,1);
  ierr = PetscUseMethod(viewer,"PetscViewerVTKGetDM_C",(PetscViewer,PetscObject*),(viewer,dm));CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PetscViewerDestroy_VTK(PetscViewer viewer)
{
  PetscViewer_VTK *vtk = (PetscViewer_VTK*)viewer->data;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscFree(vtk->filename);CHKERRQ(ierr);
  ierr = PetscFree(vtk);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)viewer,"PetscViewerFileSetName_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)viewer,"PetscViewerFileGetName_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)viewer,"PetscViewerFileSetMode_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)viewer,"PetscViewerFileGetMode_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)viewer,"PetscViewerVTKAddField_C",NULL);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)viewer,"PetscViewerVTKGetDM_C",NULL);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PetscViewerFlush_VTK(PetscViewer viewer)
{
  PetscViewer_VTK          *vtk = (PetscViewer_VTK*)viewer->data;
  PetscErrorCode           ierr;
  PetscViewerVTKObjectLink link,next;

  PetscFunctionBegin;
  if (vtk->link && (!vtk->dm || !vtk->write)) SETERRQ(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_WRONGSTATE,"No fields or no grid");
  if (vtk->write) {ierr = (*vtk->write)(vtk->dm,viewer);CHKERRQ(ierr);}
  for (link=vtk->link; link; link=next) {
    next = link->next;
    ierr = PetscObjectDestroy(&link->vec);CHKERRQ(ierr);
    ierr = PetscFree(link);CHKERRQ(ierr);
  }
  ierr       = PetscObjectDestroy(&vtk->dm);CHKERRQ(ierr);
  vtk->write = NULL;
  vtk->link  = NULL;
  PetscFunctionReturn(0);
}

PetscErrorCode  PetscViewerFileSetName_VTK(PetscViewer viewer,const char name[])
{
  PetscViewer_VTK *vtk = (PetscViewer_VTK*)viewer->data;
  PetscErrorCode  ierr;
  PetscBool       isvtk,isvts,isvtu,isvtr;
  size_t          len;

  PetscFunctionBegin;
  ierr = PetscViewerFlush(viewer);CHKERRQ(ierr);
  ierr = PetscFree(vtk->filename);CHKERRQ(ierr);
  ierr = PetscStrlen(name,&len);CHKERRQ(ierr);
  if (!len) {
    isvtk = PETSC_TRUE;
  } else {
    ierr = PetscStrcasecmp(name+len-4,".vtk",&isvtk);CHKERRQ(ierr);
    ierr = PetscStrcasecmp(name+len-4,".vts",&isvts);CHKERRQ(ierr);
    ierr = PetscStrcasecmp(name+len-4,".vtu",&isvtu);CHKERRQ(ierr);
    ierr = PetscStrcasecmp(name+len-4,".vtr",&isvtr);CHKERRQ(ierr);
  }
  if (isvtk) {
    if (viewer->format == PETSC_VIEWER_DEFAULT) viewer->format = PETSC_VIEWER_ASCII_VTK;
    if (viewer->format != PETSC_VIEWER_ASCII_VTK) SETERRQ2(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_INCOMP,"Cannot use file '%s' with format %s, should have '.vtk' extension",name,PetscViewerFormats[viewer->format]);
  } else if (isvts) {
    if (viewer->format == PETSC_VIEWER_DEFAULT) viewer->format = PETSC_VIEWER_VTK_VTS;
    if (viewer->format != PETSC_VIEWER_VTK_VTS) SETERRQ2(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_INCOMP,"Cannot use file '%s' with format %s, should have '.vts' extension",name,PetscViewerFormats[viewer->format]);
  } else if (isvtu) {
    if (viewer->format == PETSC_VIEWER_DEFAULT) viewer->format = PETSC_VIEWER_VTK_VTU;
    if (viewer->format != PETSC_VIEWER_VTK_VTU) SETERRQ2(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_INCOMP,"Cannot use file '%s' with format %s, should have '.vtu' extension",name,PetscViewerFormats[viewer->format]);
  } else if (isvtr) {
    if (viewer->format == PETSC_VIEWER_DEFAULT) viewer->format = PETSC_VIEWER_VTK_VTR;
    if (viewer->format != PETSC_VIEWER_VTK_VTR) SETERRQ2(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_INCOMP,"Cannot use file '%s' with format %s, should have '.vtr' extension",name,PetscViewerFormats[viewer->format]);
  } else SETERRQ1(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_UNKNOWN_TYPE,"File '%s' has unrecognized extension",name);
  ierr = PetscStrallocpy(len ? name : "stdout",&vtk->filename);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode  PetscViewerFileGetName_VTK(PetscViewer viewer,const char **name)
{
  PetscViewer_VTK *vtk = (PetscViewer_VTK*)viewer->data;
  PetscFunctionBegin;
  *name = vtk->filename;
  PetscFunctionReturn(0);
}

PetscErrorCode  PetscViewerFileSetMode_VTK(PetscViewer viewer,PetscFileMode type)
{
  PetscViewer_VTK *vtk = (PetscViewer_VTK*)viewer->data;

  PetscFunctionBegin;
  vtk->btype = type;
  PetscFunctionReturn(0);
}

PetscErrorCode  PetscViewerFileGetMode_VTK(PetscViewer viewer,PetscFileMode *type)
{
  PetscViewer_VTK *vtk = (PetscViewer_VTK*)viewer->data;

  PetscFunctionBegin;
  *type = vtk->btype;
  PetscFunctionReturn(0);
}

PetscErrorCode  PetscViewerVTKAddField_VTK(PetscViewer viewer,PetscObject dm,PetscErrorCode (*PetscViewerVTKWriteFunction)(PetscObject,PetscViewer),PetscInt fieldnum,PetscViewerVTKFieldType fieldtype,PetscBool checkdm,PetscObject vec)
{
  PetscViewer_VTK          *vtk = (PetscViewer_VTK*)viewer->data;
  PetscViewerVTKObjectLink link, tail = vtk->link;
  PetscErrorCode           ierr;

  PetscFunctionBegin;
  if (vtk->dm) {
    if (checkdm && dm != vtk->dm) SETERRQ(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_INCOMP,"Refusing to write a field from more than one grid to the same VTK file. Set checkdm = PETSC_FALSE to skip this check.");
  } else {
    ierr = PetscObjectReference(dm);CHKERRQ(ierr);
    vtk->dm = dm;
  }
  vtk->write  = PetscViewerVTKWriteFunction;
  ierr        = PetscNew(&link);CHKERRQ(ierr);
  link->ft    = fieldtype;
  link->vec   = vec;
  link->field = fieldnum;
  link->next  = NULL;
  /* Append to list */
  if (tail) {
    while (tail->next) tail = tail->next;
    tail->next = link;
  } else vtk->link = link;
  PetscFunctionReturn(0);
}

PetscErrorCode PetscViewerVTKGetDM_VTK(PetscViewer viewer,PetscObject *dm)
{
  PetscViewer_VTK *vtk = (PetscViewer_VTK*)viewer->data;

  PetscFunctionBegin;
  *dm = vtk->dm;
  PetscFunctionReturn(0);
}

/*MC
   PETSCVIEWERVTK - A viewer that writes to a VTK file


.seealso:  PetscViewerVTKOpen(), PetscViewerHDF5Open(), PetscViewerStringSPrintf(), PetscViewerSocketOpen(), PetscViewerDrawOpen(), PETSCVIEWERSOCKET,
           PetscViewerCreate(), PetscViewerASCIIOpen(), PetscViewerBinaryOpen(), PETSCVIEWERBINARY, PETSCVIEWERDRAW, PETSCVIEWERSTRING,
           PetscViewerMatlabOpen(), VecView(), DMView(), PetscViewerMatlabPutArray(), PETSCVIEWERASCII, PETSCVIEWERMATLAB,
           PetscViewerFileSetName(), PetscViewerFileSetMode(), PetscViewerFormat, PetscViewerType, PetscViewerSetType()

  Level: beginner
M*/

PETSC_EXTERN PetscErrorCode PetscViewerCreate_VTK(PetscViewer v)
{
  PetscViewer_VTK *vtk;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(v,&vtk);CHKERRQ(ierr);

  v->data         = (void*)vtk;
  v->ops->destroy = PetscViewerDestroy_VTK;
  v->ops->flush   = PetscViewerFlush_VTK;
  vtk->btype      = (PetscFileMode) -1;
  vtk->filename   = NULL;

  ierr = PetscObjectComposeFunction((PetscObject)v,"PetscViewerFileSetName_C",PetscViewerFileSetName_VTK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"PetscViewerFileGetName_C",PetscViewerFileGetName_VTK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"PetscViewerFileSetMode_C",PetscViewerFileSetMode_VTK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"PetscViewerFileGetMode_C",PetscViewerFileGetMode_VTK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"PetscViewerVTKAddField_C",PetscViewerVTKAddField_VTK);CHKERRQ(ierr);
  ierr = PetscObjectComposeFunction((PetscObject)v,"PetscViewerVTKGetDM_C",PetscViewerVTKGetDM_VTK);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PetscViewerVTKOpen - Opens a file for VTK output.

   Collective

   Input Parameters:
+  comm - MPI communicator
.  name - name of file
-  type - type of file
$    FILE_MODE_WRITE - create new file for binary output
$    FILE_MODE_READ - open existing file for binary input (not currently supported)
$    FILE_MODE_APPEND - open existing file for binary output (not currently supported)

   Output Parameter:
.  vtk - PetscViewer for VTK input/output to use with the specified file

   Level: beginner

   Note:
   This PetscViewer should be destroyed with PetscViewerDestroy().


.seealso: PetscViewerASCIIOpen(), PetscViewerPushFormat(), PetscViewerDestroy(),
          VecView(), MatView(), VecLoad(), MatLoad(),
          PetscFileMode, PetscViewer
@*/
PetscErrorCode PetscViewerVTKOpen(MPI_Comm comm,const char name[],PetscFileMode type,PetscViewer *vtk)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscViewerCreate(comm,vtk);CHKERRQ(ierr);
  ierr = PetscViewerSetType(*vtk,PETSCVIEWERVTK);CHKERRQ(ierr);
  ierr = PetscViewerFileSetMode(*vtk,type);CHKERRQ(ierr);
  ierr = PetscViewerFileSetName(*vtk,name);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
   PetscViewerVTKFWrite - write binary data preceded by 32-bit int length (in bytes), does not do byte swapping.

   Logically collective on PetscViewer

   Input Parameters:
+  viewer - logically collective viewer, data written from rank 0
.  fp - file pointer valid on rank 0
.  data - data pointer valid on rank 0
.  n - number of data items
-  dtype - data type

   Level: developer

   Notes:
    If PetscScalar is __float128 then the binary files are written in double precision


.seealso: DMDAVTKWriteAll(), DMComplexVTKWriteAll(), PetscViewerPushFormat(), PetscViewerVTKOpen(), PetscBinaryWrite()
@*/
PetscErrorCode PetscViewerVTKFWrite(PetscViewer viewer,FILE *fp,const void *data,PetscInt n,MPI_Datatype dtype)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;
  MPI_Datatype   vdtype=dtype;
#if defined(PETSC_USE_REAL___FLOAT128)
  double         *tmp;
  PetscInt       i;
  PetscReal      *ttmp = (PetscReal*)data;
#endif

  PetscFunctionBegin;
  if (n < 0) SETERRQ1(PetscObjectComm((PetscObject)viewer),PETSC_ERR_ARG_OUTOFRANGE,"Trying to write a negative amount of data %D",n);
  if (!n) PetscFunctionReturn(0);
  ierr = MPI_Comm_rank(PetscObjectComm((PetscObject)viewer),&rank);CHKERRQ(ierr);
  if (!rank) {
    size_t      count;
    PetscMPIInt dsize;
    PetscVTKInt bytes;

#if defined(PETSC_USE_REAL___FLOAT128)
    if (dtype == MPIU___FLOAT128) {
      ierr = PetscMalloc1(n,&tmp);CHKERRQ(ierr);
      for (i=0; i<n; i++) tmp[i] = ttmp[i];
      data  = (void*) tmp;
      vdtype = MPI_DOUBLE;
    }
#endif
    ierr  = MPI_Type_size(vdtype,&dsize);CHKERRQ(ierr);
    bytes = PetscVTKIntCast(dsize*n);

    count = fwrite(&bytes,sizeof(int),1,fp);
    if (count != 1) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_WRITE,"Error writing byte count");
    count = fwrite(data,dsize,(size_t)n,fp);
    if ((PetscInt)count != n) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_FILE_WRITE,"Wrote %D/%D array members of size %d",(PetscInt)count,n,dsize);
#if defined(PETSC_USE_REAL___FLOAT128)
    if (dtype == MPIU___FLOAT128) {
      ierr = PetscFree(tmp);CHKERRQ(ierr);
    }
#endif
  }
  PetscFunctionReturn(0);
}
