
#include <petsc/private/vecimpl.h> /*I "petscvec.h" I*/
#include "../src/vec/vec/utils/tagger/impls/andor.h"

/*@C
  VecTaggerAndGetSubs - Get the sub VecTaggers whose intersection defines the outer VecTagger

  Not collective

  Input Arguments:
. tagger - the VecTagger context

  Output Arguments:
+ nsubs - the number of sub VecTaggers
- subs - the sub VecTaggers

  Level: advanced

.seealso: VecTaggerAndSetSubs()
@*/
PetscErrorCode VecTaggerAndGetSubs(VecTagger tagger, PetscInt *nsubs, VecTagger **subs)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecTaggerGetSubs_AndOr(tagger,nsubs,subs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@C
  VecTaggerAndSetSubs - Set the sub VecTaggers whose intersection defines the outer VecTagger

  Logically collective

  Input Arguments:
+ tagger - the VecTagger context
. nsubs - the number of sub VecTaggers
- subs - the sub VecTaggers

  Level: advanced

.seealso: VecTaggerAndSetSubs()
@*/
PetscErrorCode VecTaggerAndSetSubs(VecTagger tagger, PetscInt nsubs, VecTagger *subs, PetscCopyMode mode)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecTaggerSetSubs_AndOr(tagger,nsubs,subs,mode);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode VecTaggerComputeBoxes_And(VecTagger tagger,Vec vec,PetscInt *numBoxes,VecTaggerBox **boxes)
{
  PetscInt        i, bs, nsubs, *numSubBoxes, nboxes;
  VecTaggerBox    **subBoxes;
  VecTagger       *subs;
  VecTaggerBox    *bxs = NULL;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  ierr = VecTaggerGetBlockSize(tagger,&bs);CHKERRQ(ierr);
  ierr = VecTaggerOrGetSubs(tagger,&nsubs,&subs);CHKERRQ(ierr);
  ierr = PetscMalloc2(nsubs,&numSubBoxes,nsubs,&subBoxes);CHKERRQ(ierr);
  for (i = 0; i < nsubs; i++) {
    PetscErrorCode ierr2;

    ierr2 = VecTaggerComputeBoxes(subs[i],vec,&numSubBoxes[i],&subBoxes[i]);
    if (ierr2 == PETSC_ERR_SUP) { /* no support, clean up and exit */
      PetscInt j;

      for (j = 0; j < i; j++) {
        ierr = PetscFree(subBoxes[j]);CHKERRQ(ierr);
      }
      ierr = PetscFree2(numSubBoxes,subBoxes);CHKERRQ(ierr);
      SETERRQ(PetscObjectComm((PetscObject)tagger),PETSC_ERR_SUP,"Sub tagger does not support box computation");
    } else {
      CHKERRQ(ierr2);
    }
  }
  for (i = 0, nboxes = 0; i < nsubs; i++) { /* stupid O(N^3) check to intersect boxes */
    VecTaggerBox *isect;
    PetscInt j, k, l, m, n;

    n = numSubBoxes[i];
    if (!n) {
      nboxes = 0;
      ierr = PetscFree(bxs);CHKERRQ(ierr);
      break;
    }
    if (!i) {
      ierr = PetscMalloc1(n * bs, &bxs);CHKERRQ(ierr);
      for (j = 0; j < numSubBoxes[i] * bs; j++) bxs[j] = subBoxes[i][j];
      nboxes = n;
      ierr = PetscFree(subBoxes[i]);CHKERRQ(ierr);
      continue;
    }
    ierr = PetscMalloc1(n * nboxes * bs,&isect);CHKERRQ(ierr);
    for (j = 0, l = 0; j < n; j++) {
      VecTaggerBox *subBox = &subBoxes[i][j*bs];

      for (k = 0; k < nboxes; k++) {
        PetscBool    isEmpty;
        VecTaggerBox *prevBox = &bxs[bs*k];

        ierr = VecTaggerAndOrIntersect_Private(bs,prevBox,subBox,&isect[l * bs],&isEmpty);CHKERRQ(ierr);
        if (isEmpty) continue;
        for (m = 0; m < l; m++) {
          PetscBool isSub = PETSC_FALSE;

          ierr = VecTaggerAndOrIsSubBox_Private(bs,&isect[m*bs],&isect[l*bs],&isSub);CHKERRQ(ierr);
          if (isSub) break;
          ierr = VecTaggerAndOrIsSubBox_Private(bs,&isect[l*bs],&isect[m*bs],&isSub);CHKERRQ(ierr);
          if (isSub) {
            PetscInt r;

            for (r = 0; r < bs; r++) isect[m*bs + r] = isect[l * bs + r];
            break;
          }
        }
        if (m == l) l++;
      }
    }
    ierr = PetscFree(bxs);CHKERRQ(ierr);
    bxs = isect;
    nboxes = l;
    ierr = PetscFree(subBoxes[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree2(numSubBoxes,subBoxes);CHKERRQ(ierr);
  *numBoxes = nboxes;
  *boxes = bxs;
  PetscFunctionReturn(0);
}

static PetscErrorCode VecTaggerComputeIS_And(VecTagger tagger, Vec vec, IS *is)
{
  PetscInt       nsubs, i;
  VecTagger      *subs;
  IS             isectIS;
  PetscErrorCode ierr, ierr2;

  PetscFunctionBegin;
  ierr2 = VecTaggerComputeIS_FromBoxes(tagger,vec,is);
  if (ierr2 != PETSC_ERR_SUP) {
    CHKERRQ(ierr2);
    PetscFunctionReturn(0);
  }
  ierr = VecTaggerOrGetSubs(tagger,&nsubs,&subs);CHKERRQ(ierr);
  if (!nsubs) {
    ierr = ISCreateGeneral(PetscObjectComm((PetscObject)vec),0,NULL,PETSC_OWN_POINTER,is);CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = VecTaggerComputeIS(subs[0],vec,&isectIS);CHKERRQ(ierr);
  for (i = 1; i < nsubs; i++) {
    IS subIS, newIsectIS;

    ierr = VecTaggerComputeIS(subs[i],vec,&subIS);CHKERRQ(ierr);
    ierr = ISIntersect(isectIS,subIS,&newIsectIS);CHKERRQ(ierr);
    ierr = ISDestroy(&isectIS);CHKERRQ(ierr);
    ierr = ISDestroy(&subIS);CHKERRQ(ierr);
    isectIS = newIsectIS;
  }
  *is = isectIS;
  PetscFunctionReturn(0);
}

PETSC_INTERN PetscErrorCode VecTaggerCreate_And(VecTagger tagger)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecTaggerCreate_AndOr(tagger);CHKERRQ(ierr);
  tagger->ops->computeboxes = VecTaggerComputeBoxes_And;
  tagger->ops->computeis    = VecTaggerComputeIS_And;
  PetscFunctionReturn(0);
}
