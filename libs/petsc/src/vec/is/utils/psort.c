
#include <petsc/private/petscimpl.h>
#include <petscis.h> /*I "petscis.h" I*/

/* This is the bitonic merge that works on non-power-of-2 sizes found at http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm */
static PetscErrorCode PetscParallelSortInt_Bitonic_Merge(MPI_Comm comm, PetscMPIInt tag, PetscMPIInt rankStart, PetscMPIInt rankEnd, PetscMPIInt rank, PetscMPIInt n, PetscInt keys[], PetscInt buffer[], PetscBool forward)
{
  PetscInt       diff;
  PetscInt       split, mid, partner;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  diff = rankEnd - rankStart;
  if (diff <= 0) PetscFunctionReturn(0);
  if (diff == 1) {
    if (forward) {
      ierr = PetscSortInt((PetscInt) n, keys);CHKERRQ(ierr);
    } else {
      ierr = PetscSortReverseInt((PetscInt) n, keys);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  split = 1;
  while (2 * split < diff) split *= 2;
  mid = rankStart + split;
  if (rank < mid) {
    partner = rank + split;
  } else {
    partner = rank - split;
  }
  if (partner < rankEnd) {
    PetscMPIInt i;

    ierr = MPI_Sendrecv(keys, n, MPIU_INT, partner, tag, buffer, n, MPIU_INT, partner, tag, comm, MPI_STATUS_IGNORE);CHKERRQ(ierr);
    if ((rank < partner) == (forward == PETSC_TRUE)) {
      for (i = 0; i < n; i++) {
        keys[i] = (keys[i] <= buffer[i]) ? keys[i] : buffer[i];
      }
    } else {
      for (i = 0; i < n; i++) {
        keys[i] = (keys[i] > buffer[i]) ? keys[i] : buffer[i];
      }
    }
  }
  /* divide and conquer */
  if (rank < mid) {
    ierr = PetscParallelSortInt_Bitonic_Merge(comm, tag, rankStart, mid, rank, n, keys, buffer, forward);CHKERRQ(ierr);
  } else {
    ierr = PetscParallelSortInt_Bitonic_Merge(comm, tag, mid, rankEnd, rank, n, keys, buffer, forward);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

/* This is the bitonic sort that works on non-power-of-2 sizes found at http://www.iti.fh-flensburg.de/lang/algorithmen/sortieren/bitonic/oddn.htm */
static PetscErrorCode PetscParallelSortInt_Bitonic_Recursive(MPI_Comm comm, PetscMPIInt tag, PetscMPIInt rankStart, PetscMPIInt rankEnd, PetscMPIInt rank, PetscMPIInt n, PetscInt keys[], PetscInt buffer[], PetscBool forward)
{
  PetscInt       diff;
  PetscInt       mid;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  diff = rankEnd - rankStart;
  if (diff <= 0) PetscFunctionReturn(0);
  if (diff == 1) {
    if (forward) {
      ierr = PetscSortInt(n, keys);CHKERRQ(ierr);
    } else {
      ierr = PetscSortReverseInt(n, keys);CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
  }
  mid = rankStart + diff / 2;
  /* divide and conquer */
  if (rank < mid) {
    ierr = PetscParallelSortInt_Bitonic_Recursive(comm, tag, rankStart, mid, rank, n, keys, buffer, (PetscBool) !forward);CHKERRQ(ierr);
  } else {
    ierr = PetscParallelSortInt_Bitonic_Recursive(comm, tag, mid, rankEnd, rank, n, keys, buffer, forward);CHKERRQ(ierr);
  }
  /* bitonic merge */
  ierr = PetscParallelSortInt_Bitonic_Merge(comm, tag, rankStart, rankEnd, rank, n, keys, buffer, forward);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PetscParallelSortInt_Bitonic(MPI_Comm comm, PetscInt n, PetscInt keys[])
{
  PetscMPIInt size, rank, tag, mpin;
  PetscInt       *buffer;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidIntPointer(keys, 3);
  ierr = PetscCommGetNewTag(comm, &tag);CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(comm, &rank);CHKERRQ(ierr);
  ierr = PetscMPIIntCast(n, &mpin);CHKERRQ(ierr);
  ierr = PetscMalloc1(n, &buffer);CHKERRQ(ierr);
  ierr = PetscParallelSortInt_Bitonic_Recursive(comm, tag, 0, size, rank, mpin, keys, buffer, PETSC_TRUE);CHKERRQ(ierr);
  ierr = PetscFree(buffer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PetscParallelSampleSelect(PetscLayout mapin, PetscLayout mapout, PetscInt keysin[], PetscInt *outpivots[])
{
  PetscMPIInt    size, rank;
  PetscInt       *pivots, *finalpivots, i;
  PetscInt       non_empty, my_first, count;
  PetscMPIInt    *keys_per, max_keys_per;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(mapin->comm, &size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(mapin->comm, &rank);CHKERRQ(ierr);

  /* Choose P - 1 pivots that would be ideal for the distribution on this process */
  ierr = PetscMalloc1(size - 1, &pivots);CHKERRQ(ierr);
  for (i = 0; i < size - 1; i++) {
    PetscInt index;

    if (!mapin->n) {
      /* if this rank is empty, put "infinity" in the list.  each process knows how many empty
       * processes are in the layout, so those values will be ignored from the end of the sorted
       * pivots */
      pivots[i] = PETSC_MAX_INT;
    } else {
      /* match the distribution to the desired output described by mapout by getting the index
       * that is approximately the appropriate fraction through the list */
      index = ((PetscReal) mapout->range[i + 1]) * ((PetscReal) mapin->n) / ((PetscReal) mapout->N);
      index = PetscMin(index, (mapin->n - 1));
      index = PetscMax(index, 0);
      pivots[i] = keysin[index];
    }
  }
  /* sort the pivots in parallel */
  ierr = PetscParallelSortInt_Bitonic(mapin->comm, size - 1, pivots);CHKERRQ(ierr);
#if defined(PETSC_USE_DEBUG)
  {
    PetscBool sorted;

    ierr = PetscParallelSortedInt(mapin->comm, size - 1, pivots, &sorted);CHKERRQ(ierr);
    if (!sorted) SETERRQ(mapin->comm, PETSC_ERR_PLIB, "bitonic sort failed");CHKERRQ(ierr);
  }
#endif

  /* if there are Z nonempty processes, we have (P - 1) * Z real pivots, and we want to select
   * at indices Z - 1, 2*Z - 1, ... (P - 1) * Z - 1 */
  non_empty = size;
  for (i = 0; i < size; i++) if (mapout->range[i] == mapout->range[i+1]) non_empty--;
  ierr = PetscCalloc1(size + 1, &keys_per);CHKERRQ(ierr);
  my_first = -1;
  if (non_empty) {
    for (i = 0; i < size - 1; i++) {
      size_t sample = (i + 1) * non_empty - 1;
      size_t sample_rank = sample / (size - 1);

      keys_per[sample_rank]++;
      if (my_first < 0 && (PetscMPIInt) sample_rank == rank) {
        my_first = (PetscInt) (sample - sample_rank * (size - 1));
      }
    }
  }
  for (i = 0, max_keys_per = 0; i < size; i++) max_keys_per = PetscMax(keys_per[i], max_keys_per);
  ierr = PetscMalloc1(size * max_keys_per, &finalpivots);CHKERRQ(ierr);
  /* now that we know how many pivots each process will provide, gather the selected pivots at the start of the array
   * and allgather them */
  for (i = 0; i < keys_per[rank]; i++) pivots[i] = pivots[my_first + i * non_empty];
  for (i = keys_per[rank]; i < max_keys_per; i++) pivots[i] = PETSC_MAX_INT;
  ierr = MPI_Allgather(pivots, max_keys_per, MPIU_INT, finalpivots, max_keys_per, MPIU_INT, mapin->comm);
  for (i = 0, count = 0; i < size; i++) {
    PetscInt j;

    for (j = 0; j < max_keys_per; j++) {
      if (j < keys_per[i]) {
        finalpivots[count++] = finalpivots[i * max_keys_per + j];
      }
    }
  }
  *outpivots = finalpivots;
  ierr = PetscFree(keys_per);CHKERRQ(ierr);
  ierr = PetscFree(pivots);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PetscParallelRedistribute(PetscLayout map, PetscInt n, PetscInt arrayin[], PetscInt arrayout[])
{
  PetscMPIInt    size, rank;
  PetscInt       myOffset, nextOffset;
  PetscInt       i;
  PetscMPIInt    total, filled;
  PetscMPIInt    sender, nfirst, nsecond;
  PetscMPIInt    firsttag, secondtag;
  MPI_Request    firstreqrcv;
  MPI_Request    *firstreqs;
  MPI_Request    *secondreqs;
  MPI_Status     firststatus;

  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(map->comm, &size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(map->comm, &rank);CHKERRQ(ierr);
  ierr = PetscCommGetNewTag(map->comm, &firsttag);CHKERRQ(ierr);
  ierr = PetscCommGetNewTag(map->comm, &secondtag);CHKERRQ(ierr);
  myOffset = 0;
  ierr = PetscMalloc2(size, &firstreqs, size, &secondreqs);CHKERRQ(ierr);
  ierr = MPI_Scan(&n, &nextOffset, 1, MPIU_INT, MPI_SUM, map->comm);CHKERRQ(ierr);
  myOffset = nextOffset - n;
  total = map->range[rank + 1] - map->range[rank];
  if (total > 0) {
    ierr = MPI_Irecv(arrayout, total, MPIU_INT, MPI_ANY_SOURCE, firsttag, map->comm, &firstreqrcv);CHKERRQ(ierr);
  }
  for (i = 0, nsecond = 0, nfirst = 0; i < size; i++) {
    PetscInt itotal;
    PetscInt overlap, oStart, oEnd;

    itotal = map->range[i + 1] - map->range[i];
    if (itotal <= 0) continue;
    oStart = PetscMax(myOffset, map->range[i]);
    oEnd   = PetscMin(nextOffset, map->range[i + 1]);
    overlap = oEnd - oStart;
    if (map->range[i] >= myOffset && map->range[i] < nextOffset) {
      /* send first message */
      ierr = MPI_Isend(&arrayin[map->range[i] - myOffset], overlap, MPIU_INT, i, firsttag, map->comm, &(firstreqs[nfirst++]));CHKERRQ(ierr);
    } else if (overlap > 0) {
      /* send second message */
      ierr = MPI_Isend(&arrayin[oStart - myOffset], overlap, MPIU_INT, i, secondtag, map->comm, &(secondreqs[nsecond++]));CHKERRQ(ierr);
    } else if (overlap == 0 && myOffset > map->range[i] && myOffset < map->range[i + 1]) {
      /* send empty second message */
      ierr = MPI_Isend(&arrayin[oStart - myOffset], 0, MPIU_INT, i, secondtag, map->comm, &(secondreqs[nsecond++]));CHKERRQ(ierr);
    }
  }
  filled = 0;
  sender = -1;
  if (total > 0) {
    ierr = MPI_Wait(&firstreqrcv, &firststatus);CHKERRQ(ierr);
    sender = firststatus.MPI_SOURCE;
    ierr = MPI_Get_count(&firststatus, MPIU_INT, &filled);CHKERRQ(ierr);
  }
  while (filled < total) {
    PetscMPIInt mfilled;
    MPI_Status stat;

    sender++;
    ierr = MPI_Recv(&arrayout[filled], total - filled, MPIU_INT, sender, secondtag, map->comm, &stat);CHKERRQ(ierr);
    ierr = MPI_Get_count(&stat, MPIU_INT, &mfilled);CHKERRQ(ierr);
    filled += mfilled;
  }
  ierr = MPI_Waitall(nfirst, firstreqs, MPI_STATUSES_IGNORE);CHKERRQ(ierr);
  ierr = MPI_Waitall(nsecond, secondreqs, MPI_STATUSES_IGNORE);CHKERRQ(ierr);
  ierr = PetscFree2(firstreqs, secondreqs);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode PetscParallelSortInt_Samplesort(PetscLayout mapin, PetscLayout mapout, PetscInt keysin[], PetscInt keysout[])
{
  PetscMPIInt    size, rank;
  PetscInt       *pivots = NULL, *buffer;
  PetscInt       i, j;
  PetscMPIInt    *keys_per_snd, *keys_per_rcv, *offsets_snd, *offsets_rcv, nrecv;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = MPI_Comm_size(mapin->comm, &size);CHKERRQ(ierr);
  ierr = MPI_Comm_rank(mapin->comm, &rank);CHKERRQ(ierr);
  ierr = PetscMalloc4(size, &keys_per_snd, size, &keys_per_rcv, size + 1, &offsets_snd, size + 1, &offsets_rcv);CHKERRQ(ierr);
  /* sort locally */
  ierr = PetscSortInt(mapin->n, keysin);CHKERRQ(ierr);
  /* get P - 1 pivots */
  ierr = PetscParallelSampleSelect(mapin, mapout, keysin, &pivots);CHKERRQ(ierr);
  /* determine which entries in the sorted array go to which other processes based on the pivots */
  for (i = 0, j = 0; i < size - 1; i++) {
    PetscInt prev = j;

    while ((j < mapin->n) && (keysin[j] < pivots[i])) j++;
    offsets_snd[i] = prev;
    keys_per_snd[i] = j - prev;
  }
  offsets_snd[size - 1] = j;
  keys_per_snd[size - 1] = mapin->n - j;
  offsets_snd[size] = mapin->n;
  /* get the incoming sizes */
  ierr = MPI_Alltoall(keys_per_snd, 1, MPI_INT, keys_per_rcv, 1, MPI_INT, mapin->comm);CHKERRQ(ierr);
  offsets_rcv[0] = 0;
  for (i = 0; i < size; i++) {
    offsets_rcv[i+1] = offsets_rcv[i] + keys_per_rcv[i];
  }
  nrecv = offsets_rcv[size];
  /* all to all exchange */
  ierr = PetscMalloc1(nrecv, &buffer);CHKERRQ(ierr);
  ierr = MPI_Alltoallv(keysin, keys_per_snd, offsets_snd, MPIU_INT, buffer, keys_per_rcv, offsets_rcv, MPIU_INT, mapin->comm);CHKERRQ(ierr);
  ierr = PetscFree(pivots);CHKERRQ(ierr);
  ierr = PetscFree4(keys_per_snd, keys_per_rcv, offsets_snd, offsets_rcv);CHKERRQ(ierr);

  /* local sort */
  ierr = PetscSortInt(nrecv, buffer);CHKERRQ(ierr);
#if defined(PETSC_USE_DEBUG)
  {
    PetscBool sorted;

    ierr = PetscParallelSortedInt(mapin->comm, nrecv, buffer, &sorted);CHKERRQ(ierr);
    if (!sorted) SETERRQ(mapin->comm, PETSC_ERR_PLIB, "samplesort (pre-redistribute) sort failed");CHKERRQ(ierr);
  }
#endif

  /* redistribute to the desired order */
  ierr = PetscParallelRedistribute(mapout, nrecv, buffer, keysout);CHKERRQ(ierr);
  ierr = PetscFree(buffer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*@
  PetscParallelSortInt - Globally sort a distributed array of integers

  Collective

  Input Parameters:
+ mapin - PetscLayout describing the distribution of the input keys
. mapout - PetscLayout describing the distribution of the output keys
- keysin - the pre-sorted array of integers

  Output Parameter:
. keysout - the array in which the sorted integers will be stored.  If mapin == mapout, then keysin may be equal to keysout.

  Level: developer

  Notes: This implements a distributed samplesort, which, with local array sizes n_in and n_out, global size N, and global number of processes P, does:

  - sorts locally
  - chooses pivots by sorting (in parallel) (P-1) pivot suggestions from each process using bitonic sort and allgathering a subset of (P-1) of those
  - using to the pivots to repartition the keys by all-to-all exchange
  - sorting the repartitioned keys locally (the array is now globally sorted, but does not match the mapout layout)
  - redistributing to match the mapout layout

  If keysin != keysout, then keysin will not be changed during PetscParallelSortInt.

.seealso: PetscParallelSortedInt()
@*/
PetscErrorCode PetscParallelSortInt(PetscLayout mapin, PetscLayout mapout, PetscInt keysin[], PetscInt keysout[])
{
  PetscMPIInt    size;
  PetscMPIInt    result;
  PetscInt       *keysincopy = NULL;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  PetscValidPointer(mapin, 1);
  PetscValidPointer(mapout, 2);
  ierr = MPI_Comm_compare(mapin->comm, mapout->comm, &result);CHKERRQ(ierr);
  if (result != MPI_IDENT && result != MPI_CONGRUENT) SETERRQ(mapin->comm, PETSC_ERR_ARG_NOTSAMECOMM, "layouts are not on the same communicator");
  ierr = PetscLayoutSetUp(mapin);CHKERRQ(ierr);
  ierr = PetscLayoutSetUp(mapout);CHKERRQ(ierr);
  if (mapin->n) PetscValidIntPointer(keysin, 3);
  if (mapout->n) PetscValidIntPointer(keysout, 4);
  if (mapin->N != mapout->N) SETERRQ2(mapin->comm, PETSC_ERR_ARG_SIZ, "Input and output layouts have different global sizes (%D != %D)", mapin->N, mapout->N);
  ierr = MPI_Comm_size(mapin->comm, &size);CHKERRQ(ierr);
  if (size == 1) {
    if (keysout != keysin) {
      ierr = PetscMemcpy(keysout, keysin, mapin->n * sizeof(PetscInt));CHKERRQ(ierr);
    }
    ierr = PetscSortInt(mapout->n, keysout);CHKERRQ(ierr);
    if (size == 1) PetscFunctionReturn(0);
  }
  if (keysout != keysin) {
    ierr = PetscMalloc1(mapin->n, &keysincopy);CHKERRQ(ierr);
    ierr = PetscMemcpy(keysincopy, keysin, mapin->n * sizeof(PetscInt));CHKERRQ(ierr);
    keysin = keysincopy;
  }
  ierr = PetscParallelSortInt_Samplesort(mapin, mapout, keysin, keysout);CHKERRQ(ierr);
#if defined(PETSC_USE_DEBUG)
  {
    PetscBool sorted;

    ierr = PetscParallelSortedInt(mapout->comm, mapout->n, keysout, &sorted);CHKERRQ(ierr);
    if (!sorted) SETERRQ(mapout->comm, PETSC_ERR_PLIB, "samplesort sort failed");CHKERRQ(ierr);
  }
#endif
  ierr = PetscFree(keysincopy);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
