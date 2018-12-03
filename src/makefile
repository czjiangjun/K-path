.SUFFIXES:
.SUFFIXES: .mod .o .f90
.PRECIOUS: %.f90
# .SECONDEXPANSION:

FUFFIX=.F
SUFFIX=.f90
SUCFIX=.c

OFLAG_1=-O1
OFLAG_2=-O2
OFLAG_3=-O3
OFLAG_4=

OFLAG=$(OFLAG_2)
OFLAG_IN=$(OFLAG)

LIB=lib
LLIB=-Llib -ldmy

SRCDIR=../../src
BINDIR=../../bin

include .objects
include makefile.include

# FPP=$(CPP) $*$(FUFFIX) $*$(SUFFIX) $(CPP_OPTIONS)

FPP=$(filter-out -DwNGZhalf -DNGZhalf -DwNGXhalf -DNGXhalf,$(CPP))

OBJS=$(SOURCE) $(OBJECTS)

OBJECTC = sgroup.c math_sg.c sto.c io.c lat.c pgrp.c \
pgrp_op.c  rotb.c sgrp_cub.c sgrp_hex.c sgrp_tet.c  sgrp_ort.c sgrp_trg.c \
sgrp_mon.c sgrp_trc.c

OBJSC    = sgroup.o math_sg.o sto.o io.o lat.o pgrp.o \
pgrp_op.o  rotb.o sgrp_cub.o sgrp_hex.o sgrp_tet.o  sgrp_ort.o sgrp_trg.o \
sgrp_mon.o sgrp_trc.o 

OBJS_O1=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_O1)) $(OBJECTS_O1)
OBJS_O2=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_O2)) $(OBJECTS_O2)
OBJS_O3=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_O3)) $(OBJECTS_O3)
OBJS_IN=$(filter-out $(OBJECTS_O1) $(OBJECTS_O2) $(OBJECTS_O3) $(OBJECTS_IN),$(SOURCE_IN)) $(OBJECTS_IN)

F90SRC_O1=$(OBJS_O1:.o=$(SUFFIX))
F90SRC_O2=$(OBJS_O2:.o=$(SUFFIX))
F90SRC_O3=$(OBJS_O3:.o=$(SUFFIX))
F90SRC_IN=$(OBJS_IN:.o=$(SUFFIX))

DMPI=$(findstring DMPI,$(subst DMPI_,,$(CPP)))

LIB+=$(filter-out lib,$(LIBS))
LLIB+=$(filter-out -Llib -ldmy,$(LLIBS))


ifeq ($(FFLAG_1),)
   FFLAG_1=$(FFLAG)
endif
ifeq ($(FFLAG_2),)
   FFLAG_2=$(FFLAG)
endif
ifeq ($(FFLAG_3),)
   FFLAG_3=$(FFLAG)
endif
ifeq ($(FFLAG_IN),)
   FFLAG_IN=$(FFLAG)
endif

ifeq ($(OFLAG_1),)
   OFLAG_1=$(OFLAG)
endif
ifeq ($(OFLAG_2),)
   OFLAG_2=$(OFLAG)
endif
ifeq ($(OFLAG_3),)
   OFLAG_3=$(OFLAG)
endif

ifeq ($(INCS_1),)
   INCS_1=$(INCS)
endif
ifeq ($(INCS_2),)
   INCS_2=$(INCS)
endif
ifeq ($(INCS_3),)
   INCS_3=$(INCS)
endif
ifeq ($(INCS_IN),)
   INCS_IN=$(INCS)
endif

ifndef VERSION
   ifeq ($(DMPI),DMPI)
      FPP+=-DNGZhalf
   else
      FPP+=-DNGXhalf
   endif
   EXE=kpoint_path
else
ifeq ($(VERSION),std)
   ifeq ($(DMPI),DMPI)
      FPP+=-DNGZhalf
   else
      FPP+=-DNGXhalf
   endif
   EXE=kpoint_path
else
   $(error $(VERSION) is not a valid target)
endif
endif


OBJCTS=$(filter %.o, $(OBJS) $(FFT3D))
#OBJCTS_f90=$(filter-out getshmem.o, $(OBJCTS))
F90SRC=$(OBJCTS:.o=$(SUFFIX))
F90SRC+=kpoint_path.f90

-include .depend

# export

.PHONY: all cleanall clean sources dependencies depend libs $(LIB) c.o

all: libs sources c.o
	$(MAKE) kpoint_path ; cp kpoint_path $(BINDIR)/$(EXE)
	
kpoint_path: $(OBJS) $(FFT3D) $(INC) kpoint_path.o
	rm -f kpoint_path
	$(FCL) -o kpoint_path $(OBJS) kpoint_path.o $(OBJSC) $(FFT3D) $(LLIB) $(LINK)

sources:
	rsync -u $(SRCDIR)/*.F $(SRCDIR)/*.inc $(SRCDIR)/*.c $(SRCDIR)/*.h .
	-rsync -u $(SRCDIR)/*.f .

libs: $(LIB)

$(LIB):
	rsync -ru $(SRCDIR)/$@ .
	cp makefile.include $@
	$(MAKE) -C $@
#	$(MAKE) -C $@ -j1
#
c.o: 
	$(CC_LIB) $(CFLAGS_LIB) -c $(OBJECTC)

dependencies: sources
	$(MAKE) depend

depend: $(F90SRC)
	makedepf90 $(F90SRC) > .depend

%$(SUFFIX): %$(FUFFIX)
	$(FPP)

%.o: %$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS) $(OFLAG) $(INCS) -c $*$(SUFFIX)

$(OBJS_O1): $(F90SRC_O1)
	$(FC) $(FREE) $(FFLAGS_1) $(OFLAG_1) $(INCS_1) -c $*$(SUFFIX)

$(OBJS_O2): $(F90SRC_O2)
	$(FC) $(FREE) $(FFLAGS_2) $(OFLAG_2) $(INCS_2) -c $*$(SUFFIX)

$(OBJS_O3): $(F90SRC_O3)
	$(FC) $(FREE) $(FFLAGS_3) $(OFLAG_3) $(INCS_3) -c $*$(SUFFIX)

$(OBJS_IN): $(F90SRC_IN)
	$(FC) $(FREE) $(FFLAGS_IN) $(OFLAG_IN) $(INCS_IN) -c $*$(SUFFIX)

cleanall: clean
	rm -rf $(LIB) *.F test.c

clean:
	rm -f *.o *.mod *.f90

#-----------------------------------------------------------------------
# special rules
#-----------------------------------------------------------------------

main.o : main$(SUFFIX)
	$(FC) $(FREE) $(FFLAGS) $(DEBUG) $(INCS) -c $*$(SUFFIX)

fft3dlib_f77.o: fft3dlib_f77$(SUFFIX)
	$(F77) $(FFLAGS_F77) -c $*$(SUFFIX)
