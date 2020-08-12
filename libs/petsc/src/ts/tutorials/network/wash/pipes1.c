static char help[] = "This example demonstrates the use of DMNetwork \n\\n";
/*
  Example: mpiexec -n <np> ./pipes1 -ts_max_steps 10
*/

#include "wash.h"
#include <petscdmplex.h>

/*
  WashNetworkDistribute - proc[0] distributes sequential wash object
   Input Parameters:
.  comm - MPI communicator
.  wash - wash context with all network data in proc[0]

   Output Parameter:
.  wash - wash context with nedge, nvertex and edgelist distributed
*/
PetscErrorCode WashNetworkDistribute(MPI_Comm comm,Wash wash)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank,size,tag=0;
  PetscInt       i,e,v,numEdges,numVertices,nedges,*eowners=NULL,estart,eend,*vtype=NULL,nvertices;
  PetscInt       *edgelist = wash->edgelist,*nvtx=NULL,*vtxDone=NULL;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

  numEdges    = wash->nedge;
  numVertices = wash->nvertex;

  /* (1) all processes get global and local number of edges */
  ierr = MPI_Bcast(&numEdges,1,MPIU_INT,0,comm);CHKERRQ(ierr);
  nedges = numEdges/size; /* local nedges */
  if (!rank) {
    nedges += numEdges - size*(numEdges/size);
  }
  wash->Nedge = numEdges;
  wash->nedge = nedges;
  /* ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] nedges %d, numEdges %d\n",rank,nedges,numEdges);CHKERRQ(ierr); */

  ierr = PetscCalloc3(size+1,&eowners,size,&nvtx,numVertices,&vtxDone);CHKERRQ(ierr);
  ierr = MPI_Allgather(&nedges,1,MPIU_INT,eowners+1,1,MPIU_INT,PETSC_COMM_WORLD);CHKERRQ(ierr);
  eowners[0] = 0;
  for (i=2; i<=size; i++) {
    eowners[i] += eowners[i-1];
  }

  estart = eowners[rank];
  eend   = eowners[rank+1];
  /* ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] own lists row %d - %d\n",rank,estart,eend);CHKERRQ(ierr); */

  /* (2) distribute row block edgelist to all processors */
  if (!rank) {
    vtype = wash->vtype;
    for (i=1; i<size; i++) {
      /* proc[0] sends edgelist to proc[i] */
      ierr = MPI_Send(edgelist+2*eowners[i],2*(eowners[i+1]-eowners[i]),MPIU_INT,i,tag,comm);CHKERRQ(ierr);

      /* proc[0] sends vtype to proc[i] */
      ierr = MPI_Send(vtype+2*eowners[i],2*(eowners[i+1]-eowners[i]),MPIU_INT,i,tag,comm);CHKERRQ(ierr);
    }
  } else {
    MPI_Status      status;
    ierr = PetscMalloc1(2*(eend-estart),&vtype);CHKERRQ(ierr);
    ierr = PetscMalloc1(2*(eend-estart),&edgelist);CHKERRQ(ierr);

    ierr = MPI_Recv(edgelist,2*(eend-estart),MPIU_INT,0,tag,comm,&status);CHKERRQ(ierr);
    ierr = MPI_Recv(vtype,2*(eend-estart),MPIU_INT,0,tag,comm,&status);CHKERRQ(ierr);
  }

  wash->edgelist = edgelist;

  /* (3) all processes get global and local number of vertices, without ghost vertices */
  if (!rank) {
    for (i=0; i<size; i++) {
      for (e=eowners[i]; e<eowners[i+1]; e++) {
        v = edgelist[2*e];
        if (!vtxDone[v]) {
          nvtx[i]++; vtxDone[v] = 1;
        }
        v = edgelist[2*e+1];
        if (!vtxDone[v]) {
          nvtx[i]++; vtxDone[v] = 1;
        }
      }
    }
  }
  ierr = MPI_Bcast(&numVertices,1,MPIU_INT,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = MPI_Scatter(nvtx,1,MPIU_INT,&nvertices,1,MPIU_INT,0,PETSC_COMM_WORLD);CHKERRQ(ierr);
  ierr = PetscFree3(eowners,nvtx,vtxDone);CHKERRQ(ierr);

  wash->Nvertex = numVertices;
  wash->nvertex = nvertices;
  wash->vtype   = vtype;
  PetscFunctionReturn(0);
}

PetscErrorCode WASHIFunction(TS ts,PetscReal t,Vec X,Vec Xdot,Vec F,void* ctx)
{
  PetscErrorCode ierr;
  Wash           wash=(Wash)ctx;
  DM             networkdm;
  Vec            localX,localXdot,localF, localXold;
  const PetscInt *cone;
  PetscInt       vfrom,vto,offsetfrom,offsetto,varoffset;
  PetscInt       v,vStart,vEnd,e,eStart,eEnd;
  PetscInt       nend,type;
  PetscBool      ghost;
  PetscScalar    *farr,*juncf, *pipef;
  PetscReal      dt;
  Pipe           pipe;
  PipeField      *pipex,*pipexdot,*juncx;
  Junction       junction;
  DMDALocalInfo  info;
  const PetscScalar *xarr,*xdotarr, *xoldarr;

  PetscFunctionBegin;
  localX    = wash->localX;
  localXdot = wash->localXdot;

  ierr = TSGetSolution(ts,&localXold);CHKERRQ(ierr);
  ierr = TSGetDM(ts,&networkdm);CHKERRQ(ierr);
  ierr = TSGetTimeStep(ts,&dt);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localF);CHKERRQ(ierr);

  /* Set F and localF as zero */
  ierr = VecSet(F,0.0);CHKERRQ(ierr);
  ierr = VecSet(localF,0.0);CHKERRQ(ierr);

  /* Update ghost values of locaX and locaXdot */
  ierr = DMGlobalToLocalBegin(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,X,INSERT_VALUES,localX);CHKERRQ(ierr);

  ierr = DMGlobalToLocalBegin(networkdm,Xdot,INSERT_VALUES,localXdot);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(networkdm,Xdot,INSERT_VALUES,localXdot);CHKERRQ(ierr);

  ierr = VecGetArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localXdot,&xdotarr);CHKERRQ(ierr);
  ierr = VecGetArrayRead(localXold,&xoldarr);CHKERRQ(ierr);
  ierr = VecGetArray(localF,&farr);CHKERRQ(ierr);

   /* junction->type == JUNCTION:
           juncf[0] = -qJ + sum(qin); juncf[1] = qJ - sum(qout)
       junction->type == RESERVOIR (upper stream):
           juncf[0] = -hJ + H0; juncf[1] = qJ - sum(qout)
       junction->type == VALVE (down stream):
           juncf[0] =  -qJ + sum(qin); juncf[1] = qJ
  */
  /* Vertex/junction initialization */
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  for (v=vStart; v<vEnd; v++) {
    ierr = DMNetworkIsGhostVertex(networkdm,v,&ghost);CHKERRQ(ierr);
    if (ghost) continue;

    ierr = DMNetworkGetComponent(networkdm,v,0,&type,(void**)&junction);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(networkdm,v,&varoffset);CHKERRQ(ierr);
    juncx      = (PipeField*)(xarr + varoffset);
    juncf      = (PetscScalar*)(farr + varoffset);

    juncf[0] = -juncx[0].q;
    juncf[1] =  juncx[0].q;

    if (junction->type == RESERVOIR) { /* upstream reservoir */
      juncf[0] = juncx[0].h - wash->H0;
    }
  }

  /* Edge/pipe */
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  for (e=eStart; e<eEnd; e++) {
    ierr = DMNetworkGetComponent(networkdm,e,0,&type,(void**)&pipe);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(networkdm,e,&varoffset);CHKERRQ(ierr);
    pipex    = (PipeField*)(xarr + varoffset);
    pipexdot = (PipeField*)(xdotarr + varoffset);
    pipef    = (PetscScalar*)(farr + varoffset);

    /* Get some data into the pipe structure: note, some of these operations
     * might be redundant. Will it consume too much time? */
    pipe->dt   = dt;
    pipe->xold = (PipeField*)(xoldarr + varoffset);

    /* Evaluate F over this edge/pipe: pipef[1], ...,pipef[2*nend] */
    ierr = DMDAGetLocalInfo(pipe->da,&info);CHKERRQ(ierr);
    ierr = PipeIFunctionLocal_Lax(&info,t,pipex,pipexdot,pipef,pipe);CHKERRQ(ierr);

    /* Get boundary values from connected vertices */
    ierr = DMNetworkGetConnectedVertices(networkdm,e,&cone);CHKERRQ(ierr);
    vfrom = cone[0]; /* local ordering */
    vto   = cone[1];
    ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

    /* Evaluate upstream boundary */
    ierr = DMNetworkGetComponent(networkdm,vfrom,0,&type,(void**)&junction);CHKERRQ(ierr);
    if (junction->type != JUNCTION && junction->type != RESERVOIR) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"junction type is not supported");
    juncx = (PipeField*)(xarr + offsetfrom);
    juncf = (PetscScalar*)(farr + offsetfrom);

    pipef[0] = pipex[0].h - juncx[0].h;
    juncf[1] -= pipex[0].q;

    /* Evaluate downstream boundary */
    ierr = DMNetworkGetComponent(networkdm,vto,0,&type,(void**)&junction);CHKERRQ(ierr);
    if (junction->type != JUNCTION && junction->type != VALVE) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP,"junction type is not supported");
    juncx = (PipeField*)(xarr + offsetto);
    juncf = (PetscScalar*)(farr + offsetto);
    nend  = pipe->nnodes - 1;

    pipef[2*nend + 1] = pipex[nend].h - juncx[0].h;
    juncf[0] += pipex[nend].q;
  }

  ierr = VecRestoreArrayRead(localX,&xarr);CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(localXdot,&xdotarr);CHKERRQ(ierr);
  ierr = VecRestoreArray(localF,&farr);CHKERRQ(ierr);

  ierr = DMLocalToGlobalBegin(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localF,ADD_VALUES,F);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localF);CHKERRQ(ierr);
  /*
   ierr = PetscPrintf(PETSC_COMM_WORLD("F:\n");CHKERRQ(ierr);
   ierr = VecView(F,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
   */
  PetscFunctionReturn(0);
}

PetscErrorCode WASHSetInitialSolution(DM networkdm,Vec X,Wash wash)
{
  PetscErrorCode ierr;
  PetscInt       k,nx,vkey,vfrom,vto,offsetfrom,offsetto;
  PetscInt       type,varoffset;
  PetscInt       e,eStart,eEnd;
  Vec            localX;
  PetscScalar    *xarr;
  Pipe           pipe;
  Junction       junction;
  const PetscInt *cone;
  const PetscScalar *xarray;

  PetscFunctionBegin;
  ierr = VecSet(X,0.0);CHKERRQ(ierr);
  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);
  ierr = VecGetArray(localX,&xarr);CHKERRQ(ierr);

  /* Edge */
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  for (e=eStart; e<eEnd; e++) {
    ierr = DMNetworkGetVariableOffset(networkdm,e,&varoffset);CHKERRQ(ierr);
    ierr = DMNetworkGetComponent(networkdm,e,0,&type,(void**)&pipe);CHKERRQ(ierr);

    /* set initial values for this pipe */
    ierr = PipeComputeSteadyState(pipe,wash->Q0,wash->H0);CHKERRQ(ierr);
    ierr = VecGetSize(pipe->x,&nx);CHKERRQ(ierr);

    ierr = VecGetArrayRead(pipe->x,&xarray);CHKERRQ(ierr);
    /* copy pipe->x to xarray */
    for (k=0; k<nx; k++) {
      (xarr+varoffset)[k] = xarray[k];
    }

    /* set boundary values into vfrom and vto */
    ierr = DMNetworkGetConnectedVertices(networkdm,e,&cone);CHKERRQ(ierr);
    vfrom = cone[0]; /* local ordering */
    vto   = cone[1];
    ierr = DMNetworkGetVariableOffset(networkdm,vfrom,&offsetfrom);CHKERRQ(ierr);
    ierr = DMNetworkGetVariableOffset(networkdm,vto,&offsetto);CHKERRQ(ierr);

    /* if vform is a head vertex: */
    ierr = DMNetworkGetComponent(networkdm,vfrom,0,&vkey,(void**)&junction);CHKERRQ(ierr);
    if (junction->type == RESERVOIR) {
      (xarr+offsetfrom)[1] = wash->H0; /* 1st H */
    }

    /* if vto is an end vertex: */
    ierr = DMNetworkGetComponent(networkdm,vto,0,&vkey,(void**)&junction);CHKERRQ(ierr);
    if (junction->type == VALVE) {
      (xarr+offsetto)[0] = wash->QL; /* last Q */
    }
    ierr = VecRestoreArrayRead(pipe->x,&xarray);CHKERRQ(ierr);
  }

  ierr = VecRestoreArray(localX,&xarr);CHKERRQ(ierr);
  ierr = DMLocalToGlobalBegin(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMLocalToGlobalEnd(networkdm,localX,ADD_VALUES,X);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);

#if 0
  PetscInt N;
  ierr = VecGetSize(X,&N);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"initial solution %d:\n",N);CHKERRQ(ierr);
  ierr = VecView(X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode TSDMNetworkMonitor(TS ts, PetscInt step, PetscReal t, Vec x, void *context)
{
  PetscErrorCode     ierr;
  DMNetworkMonitor   monitor;

  PetscFunctionBegin;
  monitor = (DMNetworkMonitor)context;
  ierr = DMNetworkMonitorView(monitor,x);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode PipesView(Vec X,DM networkdm,Wash wash)
{
  PetscErrorCode       ierr;
  Pipe                 pipe;
  PetscInt             key,Start,End;
  PetscMPIInt          rank;
  PetscInt             nx,nnodes,nidx,*idx1,*idx2,*idx1_h,*idx2_h,idx_start,i,k,k1,xstart,j1;
  Vec                  Xq,Xh,localX;
  IS                   is1_q,is2_q,is1_h,is2_h;
  VecScatter           ctx_q,ctx_h;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);

  /* get num of local and global total nnodes */
  nidx = wash->nnodes_loc;
  ierr = MPIU_Allreduce(&nidx,&nx,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&Xq);CHKERRQ(ierr);
  if (rank == 0) { /* all entries of Xq are in proc[0] */
    ierr = VecSetSizes(Xq,nx,PETSC_DECIDE);CHKERRQ(ierr);
  } else {
    ierr = VecSetSizes(Xq,0,PETSC_DECIDE);CHKERRQ(ierr);
  }
  ierr = VecSetFromOptions(Xq);CHKERRQ(ierr);
  ierr = VecSet(Xq,0.0);CHKERRQ(ierr);
  ierr = VecDuplicate(Xq,&Xh);CHKERRQ(ierr);

  ierr = DMGetLocalVector(networkdm,&localX);CHKERRQ(ierr);

  /* set idx1 and idx2 */
  ierr = PetscCalloc4(nidx,&idx1,nidx,&idx2,nidx,&idx1_h,nidx,&idx2_h);CHKERRQ(ierr);

  ierr = DMNetworkGetEdgeRange(networkdm,&Start, &End);CHKERRQ(ierr);

  ierr = VecGetOwnershipRange(X,&xstart,NULL);CHKERRQ(ierr);
  k1 = 0;
  j1 = 0;
  for (i = Start; i < End; i++) {
    ierr = DMNetworkGetComponent(networkdm,i,0,&key,(void**)&pipe);CHKERRQ(ierr);
    nnodes = pipe->nnodes;
    idx_start = pipe->id*nnodes;
    for (k=0; k<nnodes; k++) {
      idx1[k1] = xstart + j1*2*nnodes + 2*k;
      idx2[k1] = idx_start + k;

      idx1_h[k1] = xstart + j1*2*nnodes + 2*k + 1;
      idx2_h[k1] = idx_start + k;
      k1++;
    }
    j1++;
  }

  ierr = ISCreateGeneral(PETSC_COMM_SELF,nidx,idx1,PETSC_COPY_VALUES,&is1_q);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,nidx,idx2,PETSC_COPY_VALUES,&is2_q);CHKERRQ(ierr);
  ierr = VecScatterCreate(X,is1_q,Xq,is2_q,&ctx_q);CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx_q,X,Xq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx_q,X,Xq,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

  ierr = ISCreateGeneral(PETSC_COMM_SELF,nidx,idx1_h,PETSC_COPY_VALUES,&is1_h);CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_SELF,nidx,idx2_h,PETSC_COPY_VALUES,&is2_h);CHKERRQ(ierr);
  ierr = VecScatterCreate(X,is1_h,Xh,is2_h,&ctx_h);CHKERRQ(ierr);
  ierr = VecScatterBegin(ctx_h,X,Xh,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);
  ierr = VecScatterEnd(ctx_h,X,Xh,INSERT_VALUES,SCATTER_FORWARD);CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Xq: \n");CHKERRQ(ierr);
  ierr = VecView(Xq,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Xh: \n");CHKERRQ(ierr);
  ierr = VecView(Xh,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  ierr = VecScatterDestroy(&ctx_q);CHKERRQ(ierr);
  ierr = PetscFree4(idx1,idx2,idx1_h,idx2_h);CHKERRQ(ierr);
  ierr = ISDestroy(&is1_q);CHKERRQ(ierr);
  ierr = ISDestroy(&is2_q);CHKERRQ(ierr);

  ierr = VecScatterDestroy(&ctx_h);CHKERRQ(ierr);
  ierr = ISDestroy(&is1_h);CHKERRQ(ierr);
  ierr = ISDestroy(&is2_h);CHKERRQ(ierr);

  ierr = VecDestroy(&Xq);CHKERRQ(ierr);
  ierr = VecDestroy(&Xh);CHKERRQ(ierr);
  ierr = DMRestoreLocalVector(networkdm,&localX);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode WashNetworkCleanUp(Wash wash)
{
  PetscErrorCode ierr;
  PetscMPIInt    rank;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(wash->comm,&rank);CHKERRQ(ierr);
  ierr = PetscFree(wash->edgelist);CHKERRQ(ierr);
  ierr = PetscFree(wash->vtype);CHKERRQ(ierr);
  if (!rank) {
    ierr = PetscFree2(wash->junction,wash->pipe);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

PetscErrorCode WashNetworkCreate(MPI_Comm comm,PetscInt pipesCase,Wash *wash_ptr)
{
  PetscErrorCode ierr;
  PetscInt       npipes;
  PetscMPIInt    rank;
  Wash           wash=NULL;
  PetscInt       i,numVertices,numEdges,*vtype;
  PetscInt       *edgelist;
  Junction       junctions=NULL;
  Pipe           pipes=NULL;
  PetscBool      washdist=PETSC_TRUE;

  PetscFunctionBegin;
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);

  ierr = PetscCalloc1(1,&wash);CHKERRQ(ierr);
  wash->comm = comm;
  *wash_ptr  = wash;
  wash->Q0   = 0.477432; /* RESERVOIR */
  wash->H0   = 150.0;
  wash->HL   = 143.488;  /* VALVE */
  wash->QL   = 0.0;
  wash->nnodes_loc = 0;

  numVertices = 0;
  numEdges    = 0;
  edgelist    = NULL;

  /* proc[0] creates a sequential wash and edgelist */
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Setup pipesCase %D\n",pipesCase);CHKERRQ(ierr);

  /* Set global number of pipes, edges, and junctions */
  /*-------------------------------------------------*/
  switch (pipesCase) {
  case 0:
    /* pipeCase 0: */
    /* =================================================
    (RESERVOIR) v0 --E0--> v1--E1--> v2 --E2-->v3 (VALVE)
    ====================================================  */
    npipes = 3;
    ierr = PetscOptionsGetInt(NULL,NULL, "-npipes", &npipes, NULL);CHKERRQ(ierr);
    wash->nedge   = npipes;
    wash->nvertex = npipes + 1;

    /* Set local edges and vertices -- proc[0] sets entire network, then distributes */
    numVertices = 0;
    numEdges    = 0;
    edgelist    = NULL;
    if (!rank) {
      numVertices = wash->nvertex;
      numEdges    = wash->nedge;

      ierr = PetscCalloc1(2*numEdges,&edgelist);CHKERRQ(ierr);
      for (i=0; i<numEdges; i++) {
        edgelist[2*i] = i; edgelist[2*i+1] = i+1;
      }

      /* Add network components */
      /*------------------------*/
      ierr = PetscCalloc2(numVertices,&junctions,numEdges,&pipes);CHKERRQ(ierr);

      /* vertex */
      for (i=0; i<numVertices; i++) {
        junctions[i].id = i;
        junctions[i].type = JUNCTION;
      }

      junctions[0].type             = RESERVOIR;
      junctions[numVertices-1].type = VALVE;
    }
    break;
  case 1:
    /* pipeCase 1: */
    /* ==========================
                v2 (VALVE)
                ^
                |
               E2
                |
    v0 --E0--> v3--E1--> v1
  (RESERVOIR)            (RESERVOIR)
    =============================  */
    npipes = 3;
    wash->nedge   = npipes;
    wash->nvertex = npipes + 1;

    /* Set local edges and vertices -- proc[0] sets entire network, then distributes */
    if (!rank) {
      numVertices = wash->nvertex;
      numEdges    = wash->nedge;

      ierr = PetscCalloc1(2*numEdges,&edgelist);CHKERRQ(ierr);
      edgelist[0] = 0; edgelist[1] = 3;  /* edge[0] */
      edgelist[2] = 3; edgelist[3] = 1;  /* edge[1] */
      edgelist[4] = 3; edgelist[5] = 2;  /* edge[2] */

      /* Add network components */
      /*------------------------*/
      ierr = PetscCalloc2(numVertices,&junctions,numEdges,&pipes);CHKERRQ(ierr);
      /* vertex */
      for (i=0; i<numVertices; i++) {
        junctions[i].id   = i;
        junctions[i].type = JUNCTION;
      }

      junctions[0].type = RESERVOIR;
      junctions[1].type = VALVE;
      junctions[2].type = VALVE;
    }
    break;
  case 2:
    /* pipeCase 2: */
    /* ==========================
    (RESERVOIR)  v2--> E2
                       |
            v0 --E0--> v3--E1--> v1
    (RESERVOIR)               (VALVE)
    =============================  */

    /* Set application parameters -- to be used in function evalutions */
    npipes = 3;
    wash->nedge   = npipes;
    wash->nvertex = npipes + 1;

    /* Set local edges and vertices -- proc[0] sets entire network, then distributes */
    if (!rank) {
      numVertices = wash->nvertex;
      numEdges    = wash->nedge;

      ierr = PetscCalloc1(2*numEdges,&edgelist);CHKERRQ(ierr);
      edgelist[0] = 0; edgelist[1] = 3;  /* edge[0] */
      edgelist[2] = 3; edgelist[3] = 1;  /* edge[1] */
      edgelist[4] = 2; edgelist[5] = 3;  /* edge[2] */

      /* Add network components */
      /*------------------------*/
      ierr = PetscCalloc2(numVertices,&junctions,numEdges,&pipes);CHKERRQ(ierr);
      /* vertex */
      for (i=0; i<numVertices; i++) {
        junctions[i].id = i;
        junctions[i].type = JUNCTION;
      }

      junctions[0].type = RESERVOIR;
      junctions[1].type = VALVE;
      junctions[2].type = RESERVOIR;
    }
    break;
  default:
    SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"not done yet");
  }

  /* set edge global id */
  for (i=0; i<numEdges; i++) pipes[i].id = i;

  if (!rank) { /* set vtype for proc[0] */
    PetscInt v;
    ierr = PetscMalloc1(2*numEdges,&vtype);CHKERRQ(ierr);
    for (i=0; i<2*numEdges; i++) {
      v        = edgelist[i];
      vtype[i] = junctions[v].type;
    }
    wash->vtype = vtype;
  }

  *wash_ptr      = wash;
  wash->nedge    = numEdges;
  wash->nvertex  = numVertices;
  wash->edgelist = edgelist;
  wash->junction = junctions;
  wash->pipe     = pipes;

  /* Distribute edgelist to other processors */
  ierr = PetscOptionsGetBool(NULL,NULL,"-wash_distribute",&washdist,NULL);CHKERRQ(ierr);
  if (washdist) {
    /*
     ierr = PetscPrintf(PETSC_COMM_WORLD," Distribute sequential wash ...\n");CHKERRQ(ierr);
     */
    ierr = WashNetworkDistribute(comm,wash);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* ------------------------------------------------------- */
int main(int argc,char ** argv)
{
  PetscErrorCode    ierr;
  Wash              wash;
  Junction          junctions,junction;
  Pipe              pipe,pipes;
  PetscInt          KeyPipe,KeyJunction;
  PetscInt          *edgelist = NULL,*edgelists[1],*vtype = NULL;
  PetscInt          i,e,v,eStart,eEnd,vStart,vEnd,key;
  PetscInt          vkey,type;
  const PetscInt    *cone;
  DM                networkdm;
  PetscMPIInt       size,rank;
  PetscReal         ftime;
  Vec               X;
  TS                ts;
  PetscInt          steps=1;
  TSConvergedReason reason;
  PetscBool         viewpipes,monipipes=PETSC_FALSE,userJac=PETSC_TRUE,viewdm=PETSC_FALSE,viewX=PETSC_FALSE;
  PetscBool         test=PETSC_FALSE;
  PetscInt          pipesCase=0;
  DMNetworkMonitor  monitor;
  MPI_Comm          comm;

  PetscInt          nedges,nvertices; /* local num of edges and vertices */
  PetscInt          nnodes = 6,nv,ne;
  const PetscInt    *vtx,*edge;

  ierr = PetscInitialize(&argc,&argv,"pOption",help);if (ierr) return ierr;

  /* Read runtime options */
  ierr = PetscOptionsGetInt(NULL,NULL, "-case", &pipesCase, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-user_Jac",&userJac,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-pipe_monitor",&monipipes,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-viewdm",&viewdm,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-viewX",&viewX,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL, "-npipenodes", &nnodes, NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetBool(NULL,NULL,"-test",&test,NULL);CHKERRQ(ierr);

  /* Create networkdm */
  /*------------------*/
  ierr = DMNetworkCreate(PETSC_COMM_WORLD,&networkdm);CHKERRQ(ierr);
  ierr = PetscObjectGetComm((PetscObject)networkdm,&comm);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm,&size);CHKERRQ(ierr);

  if (size == 1 && monipipes) {
    ierr = DMNetworkMonitorCreate(networkdm,&monitor);CHKERRQ(ierr);
  }

  /* Register the components in the network */
  ierr = DMNetworkRegisterComponent(networkdm,"junctionstruct",sizeof(struct _p_Junction),&KeyJunction);CHKERRQ(ierr);
  ierr = DMNetworkRegisterComponent(networkdm,"pipestruct",sizeof(struct _p_Pipe),&KeyPipe);CHKERRQ(ierr);

  /* Create a distributed wash network (user-specific) */
  ierr = WashNetworkCreate(comm,pipesCase,&wash);CHKERRQ(ierr);
  nedges      = wash->nedge;
  nvertices   = wash->nvertex; /* local num of vertices, excluding ghosts */
  edgelist    = wash->edgelist;
  vtype       = wash->vtype;
  junctions   = wash->junction;
  pipes       = wash->pipe;

#if 0
  ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] after WashNetworkDistribute...ne %d, nv %d\n",rank,nedges,nvertices);CHKERRQ(ierr);
  if (rank == 0) {
    for (e=0; e<nedges; e++) ierr = PetscPrintf(PETSC_COMM_SELF," edge %d --> %d\n",edgelist[2*e],edgelist[2*e+1]);CHKERRQ(ierr);
  }
#endif

  /* Set up the network layout */
  ierr = DMNetworkSetSizes(networkdm,1,&nvertices,&nedges,0,NULL);CHKERRQ(ierr);

  /* Add local edge connectivity */
  edgelists[0] = edgelist;
  ierr = DMNetworkSetEdgeList(networkdm,edgelists,NULL);CHKERRQ(ierr);
  ierr = DMNetworkLayoutSetUp(networkdm);CHKERRQ(ierr);

  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
  /* ierr = PetscPrintf(PETSC_COMM_SELF,"[%d] eStart/End: %d - %d; vStart/End: %d - %d\n",rank,eStart,eEnd,vStart,vEnd);CHKERRQ(ierr); */

  /* Test DMNetworkGetSubnetworkInfo() */
  if (test) {
    ierr = DMNetworkGetSubnetworkInfo(networkdm,0,&nv,&ne,&vtx,&edge);CHKERRQ(ierr);
    if (ne != eEnd - eStart || nv != vEnd - vStart) SETERRQ2(PetscObjectComm((PetscObject)networkdm),PETSC_ERR_ARG_WRONG,"ne %D or nv %D is incorrect",ne,nv);
  }

  if (rank) { /* junctions[] and pipes[] for proc[0] are allocated in WashNetworkCreate() */
    /* vEnd - vStart = nvertices + num of ghost vertices! */
    ierr = PetscCalloc2(vEnd - vStart,&junctions,nedges,&pipes);CHKERRQ(ierr);
  }

  /* Add Pipe component to all local edges */
  for (e = eStart; e < eEnd; e++) {
    if (test) {
      if (e != edge[e]) SETERRQ2(PetscObjectComm((PetscObject)networkdm),PETSC_ERR_ARG_WRONG,"e %D != edge %D from DMNetworkGetSubnetworkInfo()",e,edge[e]);
    }

    pipes[e-eStart].nnodes = nnodes;
    ierr = DMNetworkAddComponent(networkdm,e,KeyPipe,&pipes[e-eStart]);CHKERRQ(ierr);

    /* Add number of variables to each edge */
    ierr = DMNetworkAddNumVariables(networkdm,e,2*pipes[e-eStart].nnodes);CHKERRQ(ierr);

    if (size == 1 && monipipes) { /* Add monitor -- show Q_{pipes[e-eStart].id}? */
      pipes[e-eStart].length = 600.0;
      ierr = DMNetworkMonitorAdd(monitor, "Pipe Q", e, pipes[e-eStart].nnodes, 0, 2, 0.0,pipes[e-eStart].length, -0.8, 0.8, PETSC_TRUE);CHKERRQ(ierr);
      ierr = DMNetworkMonitorAdd(monitor, "Pipe H", e, pipes[e-eStart].nnodes, 1, 2, 0.0,pipes[e-eStart].length, -400.0, 800.0, PETSC_TRUE);CHKERRQ(ierr);
    }
  }

  /* Add Junction component to all local vertices, including ghost vertices! */
  for (v = vStart; v < vEnd; v++) {
    if (test) {
      if (v != vtx[v-vStart]) SETERRQ2(PetscObjectComm((PetscObject)networkdm),PETSC_ERR_ARG_WRONG,"v %D != vtx %D from DMNetworkGetSubnetworkInfo()",v,vtx[v-vStart]);
    }

    ierr = DMNetworkAddComponent(networkdm,v,KeyJunction,&junctions[v-vStart]);CHKERRQ(ierr);

    /* Add number of variables to vertex */
    ierr = DMNetworkAddNumVariables(networkdm,v,2);CHKERRQ(ierr);
  }

  if (size > 1) {  /* must be called before DMSetUp()???. Other partitioners do not work yet??? -- cause crash in proc[0]! */
    DM               plexdm;
    PetscPartitioner part;
    ierr = DMNetworkGetPlex(networkdm,&plexdm);CHKERRQ(ierr);
    ierr = DMPlexGetPartitioner(plexdm, &part);CHKERRQ(ierr);
    ierr = PetscPartitionerSetType(part,PETSCPARTITIONERSIMPLE);CHKERRQ(ierr);
  }

  /* Set up DM for use */
  ierr = DMSetUp(networkdm);CHKERRQ(ierr);
  if (viewdm) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\nAfter DMSetUp, DMView:\n");CHKERRQ(ierr);
    ierr = DMView(networkdm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  /* Set user physical parameters to the components */
  for (e = eStart; e < eEnd; e++) {
    ierr = DMNetworkGetConnectedVertices(networkdm,e,&cone);CHKERRQ(ierr);
    /* vfrom */
    ierr = DMNetworkGetComponent(networkdm,cone[0],0,&vkey,(void**)&junction);CHKERRQ(ierr);
    junction->type = (VertexType)vtype[2*e];

    /* vto */
    ierr = DMNetworkGetComponent(networkdm,cone[1],0,&vkey,(void**)&junction);CHKERRQ(ierr);
    junction->type = (VertexType)vtype[2*e+1];
  }

  ierr = WashNetworkCleanUp(wash);CHKERRQ(ierr);

  /* Network partitioning and distribution of data */
  ierr = DMNetworkDistribute(&networkdm,0);CHKERRQ(ierr);
  if (viewdm) {
    PetscPrintf(PETSC_COMM_WORLD,"\nAfter DMNetworkDistribute, DMView:\n");CHKERRQ(ierr);
    ierr = DMView(networkdm,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  /* create vectors */
  ierr = DMCreateGlobalVector(networkdm,&X);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(networkdm,&wash->localX);CHKERRQ(ierr);
  ierr = DMCreateLocalVector(networkdm,&wash->localXdot);CHKERRQ(ierr);

  /* PipeSetUp -- each process only sets its own pipes */
  /*---------------------------------------------------*/
  ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);

  userJac = PETSC_TRUE;
  ierr = DMNetworkHasJacobian(networkdm,userJac,userJac);CHKERRQ(ierr);
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
  for (e=eStart; e<eEnd; e++) { /* each edge has only one component, pipe */
    ierr = DMNetworkGetComponent(networkdm,e,0,&type,(void**)&pipe);CHKERRQ(ierr);

    wash->nnodes_loc += pipe->nnodes; /* local total num of nodes, will be used by PipesView() */
    ierr = PipeSetParameters(pipe,
                             600.0,          /* length */
                             0.5,            /* diameter */
                             1200.0,         /* a */
                             0.018);CHKERRQ(ierr);    /* friction */
    ierr = PipeSetUp(pipe);CHKERRQ(ierr);

    if (userJac) {
      /* Create Jacobian matrix structures for a Pipe */
      Mat            *J;
      ierr = PipeCreateJacobian(pipe,NULL,&J);CHKERRQ(ierr);
      ierr = DMNetworkEdgeSetMatrix(networkdm,e,J);CHKERRQ(ierr);
    }
  }

  if (userJac) {
    ierr = DMNetworkGetVertexRange(networkdm,&vStart,&vEnd);CHKERRQ(ierr);
    for (v=vStart; v<vEnd; v++) {
      Mat            *J;
      ierr = JunctionCreateJacobian(networkdm,v,NULL,&J);CHKERRQ(ierr);
      ierr = DMNetworkVertexSetMatrix(networkdm,v,J);CHKERRQ(ierr);

      ierr = DMNetworkGetComponent(networkdm,v,0,&vkey,(void**)&junction);CHKERRQ(ierr);
      junction->jacobian = J;
    }
  }

  /* Test DMNetworkGetSubnetworkInfo() */
  if (test) {
    ierr = DMNetworkGetEdgeRange(networkdm,&eStart,&eEnd);CHKERRQ(ierr);
    ierr = DMNetworkGetSubnetworkInfo(networkdm,0,&nv,&ne,&vtx,&edge);CHKERRQ(ierr);
    if (ne != eEnd - eStart || nv != vEnd - vStart) SETERRQ2(PetscObjectComm((PetscObject)networkdm),PETSC_ERR_ARG_WRONG,"ne %D or nv %D is incorrect",ne,nv);

    for (e = eStart; e < eEnd; e++) {
      if (e != edge[e]) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"e %D != edge %D from DMNetworkGetSubnetworkInfo()",e,edge[e]);
    }
    for (v = vStart; v < vEnd; v++) {
      if (v != vtx[v-vStart]) SETERRQ2(PetscObjectComm((PetscObject)networkdm),PETSC_ERR_ARG_WRONG,"v %D != vtx %D from DMNetworkGetSubnetworkInfo()",v,vtx[v-vStart]);
    }
  }

  /* Setup solver                                           */
  /*--------------------------------------------------------*/
  ierr = TSCreate(PETSC_COMM_WORLD,&ts);CHKERRQ(ierr);

  ierr = TSSetDM(ts,(DM)networkdm);CHKERRQ(ierr);
  ierr = TSSetIFunction(ts,NULL,WASHIFunction,wash);CHKERRQ(ierr);

  ierr = TSSetMaxSteps(ts,steps);CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_STEPOVER);CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,0.1);CHKERRQ(ierr);
  ierr = TSSetType(ts,TSBEULER);CHKERRQ(ierr);
  if (size == 1 && monipipes) {
    ierr = TSMonitorSet(ts, TSDMNetworkMonitor, monitor, NULL);CHKERRQ(ierr);
  }
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

  ierr = WASHSetInitialSolution(networkdm,X,wash);CHKERRQ(ierr);

  ierr = TSSolve(ts,X);CHKERRQ(ierr);

  ierr = TSGetSolveTime(ts,&ftime);CHKERRQ(ierr);
  ierr = TSGetStepNumber(ts,&steps);CHKERRQ(ierr);
  ierr = TSGetConvergedReason(ts,&reason);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"%s at time %g after %D steps\n",TSConvergedReasons[reason],(double)ftime,steps);CHKERRQ(ierr);
  if (viewX) {
    ierr = VecView(X,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  }

  viewpipes = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL, "-Jac_view", &viewpipes,NULL);CHKERRQ(ierr);
  if (viewpipes) {
    SNES snes;
    Mat  Jac;
    ierr = TSGetSNES(ts,&snes);CHKERRQ(ierr);
    ierr = SNESGetJacobian(snes,&Jac,NULL,NULL,NULL);CHKERRQ(ierr);
    ierr = MatView(Jac,PETSC_VIEWER_DRAW_WORLD);CHKERRQ(ierr);
  }

  /* View solution q and h */
  /* --------------------- */
  viewpipes = PETSC_FALSE;
  ierr = PetscOptionsGetBool(NULL,NULL, "-pipe_view", &viewpipes,NULL);CHKERRQ(ierr);
  if (viewpipes) {
    ierr = PipesView(X,networkdm,wash);CHKERRQ(ierr);
  }

  /* Free spaces */
  /* ----------- */
  ierr = TSDestroy(&ts);CHKERRQ(ierr);
  ierr = VecDestroy(&X);CHKERRQ(ierr);
  ierr = VecDestroy(&wash->localX);CHKERRQ(ierr);
  ierr = VecDestroy(&wash->localXdot);CHKERRQ(ierr);

  /* Destroy objects from each pipe that are created in PipeSetUp() */
  ierr = DMNetworkGetEdgeRange(networkdm,&eStart, &eEnd);CHKERRQ(ierr);
  for (i = eStart; i < eEnd; i++) {
    ierr = DMNetworkGetComponent(networkdm,i,0,&key,(void**)&pipe);CHKERRQ(ierr);
    ierr = PipeDestroy(&pipe);CHKERRQ(ierr);
  }
  if (userJac) {
    for (v=vStart; v<vEnd; v++) {
      ierr = DMNetworkGetComponent(networkdm,v,0,&vkey,(void**)&junction);CHKERRQ(ierr);
      ierr = JunctionDestroyJacobian(networkdm,v,junction);CHKERRQ(ierr);
    }
  }

  if (size == 1 && monipipes) {
    ierr = DMNetworkMonitorDestroy(&monitor);CHKERRQ(ierr);
  }
  ierr = DMDestroy(&networkdm);CHKERRQ(ierr);
  ierr = PetscFree(wash);CHKERRQ(ierr);

  if (rank) {
    ierr = PetscFree2(junctions,pipes);CHKERRQ(ierr);
  }
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

   build:
     depends: pipeInterface.c pipeImpls.c

   test:
      args: -ts_monitor -case 1 -ts_max_steps 1 -options_left no -viewX -test
      localrunfiles: pOption
      output_file: output/pipes1_1.out

   test:
      suffix: 2
      nsize: 2
      requires: mumps
      args: -ts_monitor -case 1 -ts_max_steps 1 -petscpartitioner_type simple -options_left no -viewX -test
      localrunfiles: pOption
      output_file: output/pipes1_2.out

   test:
      suffix: 3
      nsize: 2
      requires: mumps
      args: -ts_monitor -case 0 -ts_max_steps 1 -petscpartitioner_type simple -options_left no -viewX -test
      localrunfiles: pOption
      output_file: output/pipes1_3.out

   test:
      suffix: 4
      args: -ts_monitor -case 2 -ts_max_steps 1 -options_left no -viewX -test
      localrunfiles: pOption
      output_file: output/pipes1_4.out

   test:
      suffix: 5
      nsize: 3
      requires: mumps
      args: -ts_monitor -case 2 -ts_max_steps 10 -petscpartitioner_type simple -options_left no -viewX -test
      localrunfiles: pOption
      output_file: output/pipes1_5.out

   test:
      suffix: 6
      nsize: 2
      requires: mumps
      args: -ts_monitor -case 1 -ts_max_steps 1 -petscpartitioner_type simple -options_left no -wash_distribute 0 -viewX -test
      localrunfiles: pOption
      output_file: output/pipes1_6.out

   test:
      suffix: 7
      nsize: 2
      requires: mumps
      args: -ts_monitor -case 2 -ts_max_steps 1 -petscpartitioner_type simple -options_left no -wash_distribute 0 -viewX -test
      localrunfiles: pOption
      output_file: output/pipes1_7.out

TEST*/
