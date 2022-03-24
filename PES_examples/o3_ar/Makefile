.SUFFIXES:
.SUFFIXES: .f .f90 .o

FL = ftn
CMPLFLG = -c -traceback

LNKFLG = -o $(TARGET)
TARGET = O3-Ar_PES.x

OBJS = MODULES_rigid3D.o POTEN_rigid3D.o SUBROUTINES_F77.o example.o

$(TARGET) : $(OBJS)
	$(FL) $(OBJS) $(LNKFLG)
$(OBJS) :
.f90.o:
	$(FL) $(CMPLFLG) $<

.f.o:
	$(FL) $(CMPLFLG) $<

clean:
	rm *.mod *.o O3-Ar_PES.x
  