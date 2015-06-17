
F90 = gfortran
F90OPTS = -cpp -g -ffixed-line-length-none -ffree-line-length-none -fopenmp -O2
FLAGS=-lblas -llapack

PROJECTION=iomod.o \
	projection.o

MKINFO=iomod.o \
	rddatamod.o \
	eckartmod.o \
	mkinfo.o

DISPLACE=iomod.o \
	rddatamod.o \
	displace.o

VIBCP=iomod.o \
	rddatamod.o \
	vibcp.o

ALL=iomod.o \
	rddatamod.o \
	eckartmod.o \
	projection.o \
	mkinfo.o \
	displace.o \
	vibcp.o

projection: $(PROJECTION)
	$(F90) $(F90OPTS) $(PROJECTION) -o projection $(FLAGS)
	rm -f *.o *~ *.mod

mkinfo: $(MKINFO)
	$(F90) $(F90OPTS) $(MKINFO) -o mkinfo $(FLAGS)
	rm -f *.o *~ *.mod

displace: $(DISPLACE)
	$(F90) $(F90OPTS) $(DISPLACE) -o displace $(FLAGS)
	rm -f *.o *~ *.mod

vibcp: $(VIBCP)
	$(F90) $(F90OPTS) $(VIBCP) -o vibcp $(FLAGS)
	rm -f *.o *~ *.mod

all: $(ALL)
	$(F90) $(F90OPTS) $(PROJECTION) -o projection $(FLAGS)
	$(F90) $(F90OPTS) $(MKINFO) -o mkinfo $(FLAGS)
	$(F90) $(F90OPTS) $(DISPLACE) -o displace $(FLAGS)
	$(F90) $(F90OPTS) $(VIBCP) -o vibcp $(FLAGS)
	rm -f *.o *~ *.mod

%.o: %.f90
	$(F90) -c $(F90OPTS) $<

%.o: %.F90
	$(F90) -c $(F90OPTS) $<

