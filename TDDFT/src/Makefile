# Makefile for ce-tddft

QE_SOURCE=/Users/twh/Documents/espresso/
include $(QE_SOURCE)/make.sys
DFLAGS += 

# location of needed modules
MODFLAGS = $(MOD_FLAG)$(QE_SOURCE)/iotk/src $(MOD_FLAG)$(QE_SOURCE)/Modules \
           $(MOD_FLAG)$(QE_SOURCE)/PW/src $(MOD_FLAG).

TDDFT_OBJS = tddft_version.o tddft_module.o ruku.o  cn.o tddft_routines.o \
	     tddft_setup.o tddft_ch_psi_all.o tddft_cgsolver.o \
	     update_ham.o apply_efield.o molecule_optical_absorption.o \
	     molecule_operators.o ehrenfest_routines.o apply_efield_t.o \
	     add_efield_t.o fparser.o sum_band_tddft.o crank_nicolson.o \
           runge_kutta.o sinv_hpsi.o hpsi_deltat.o vhpsi_deltat.o

LIBOBJS = $(QE_SOURCE)/flib/ptools.a $(QE_SOURCE)/flib/flib.a $(QE_SOURCE)/clib/clib.a \
 	  $(QE_SOURCE)/iotk/src/libiotk.a
ifeq ($(wildcard $(QE_SOURCE)/Environ),)
  QEMODS = $(QE_SOURCE)/Modules/libqemod.a
else
  QEMODS = $(QE_SOURCE)/Environ/src/libenviron.a $(QE_SOURCE)/Modules/libqemod.a
endif
PWOBJS  = $(QE_SOURCE)/PW/src/libpw.a

all: version tddft.x fdetddft.x

version:
	./update_tddft_version

tddft_version.f90: version

tddft.x: $(TDDFT_OBJS) tddft_main.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ $(TDDFT_OBJS) tddft_main.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	(cd ../bin; ln -sf ../src/tddft.x .)
	(test -d ../../bin || exit 1; cd ../../bin; ln -sf ../TDDFT/src/tddft.x .)

fdetddft.x: $(TDDFT_OBJS) tddft_main_fde.o $(PWOBJS) $(QEMODS) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ $(TDDFT_OBJS) tddft_main_fde.o $(PWOBJS) $(QEMODS) $(LIBOBJS) $(LIBS)
	(cd ../bin; ln -sf ../src/fdetddft.x .)
	(test -d ../../bin || exit 1; cd ../../bin; ln -sf ../TDDFT/src/fdetddft.x .)

clean:
	-/bin/rm -f tddft.x fdetddft.x *.o *.F90 *__genmod.f90 *.d *.mod *.i *.L

distclean: clean
	-/bin/rm -f Makefile make.depend

include make.depend

# DO NOT DELETE

