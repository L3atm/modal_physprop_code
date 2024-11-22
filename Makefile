#MAKEFILE for modal optics code

#choose options based on compilers
ifeq ($(COMP) , lf95)
#For Lahey compiler
FC = lf95
LDFLAGS = 	
LDLIBS = 
FCFLAGS = 
ifeq ($(DEBUG) , y)
FCFLAGS = --chk a,e,s,u --trace --trap -g --O0
endif
NETCDF_INC = /share/apps/netcdf/4.1.3/lahey/8.1/include
NETCDF_LIB = /share/apps/netcdf/4.1.3/lahey/8.1/lib
endif

ifeq ($(COMP) , ifort)
#For Intel compiler	
FC = ifort	
LDFLAGS = 	
LDLIBS = 
FCFLAGS = -fpp
ifeq ($(DEBUG) , y)
FCFLAGS = -ftz -convert big_endian -fp-model precise  -check all -check noarg_temp_created -fpe0 -g -traceback 
endif
#hbrown, changed NETCDF paths
NETCDF_INC = /gpfs/u/apps/opt/netcdf-mpi-tdsf/4.3.2/intel/12.1.5/include
NETCDF_LIB = /gpfs/u/apps/opt/netcdf-mpi-tdsf/4.3.2/intel/12.1.5/lib
endif


CPP = cpp
CPPFLAGS = $(mam_modes) -DMODAL_AERO
EXE = optics.exe

#
# Suffix rules...
#
.SUFFIXES : .F .f .F90 .f90 .o

.F90.o:
	$(RM) $*.mod ; $(CPP) $(CPPFLAGS) $*.F90 $*.f90 ;  $(FC) -c $(FCFLAGS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdf -lnetcdff  $*.f90

.F.o:
	$(RM) $*.mod ; $(FC) -c $(FCFLAGS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdf -lnetcdff   $*.F


#
# Object lists...
#

OBJ =                           \
	shr_kind_mod.o          \
	shr_const_mod.o         \
	physconst.o             \
	gen_modal_optics_file.o   \
	modal_aero_data.o       \
	radconstants.o          \
	module_optics_support.o \
	miesubs.o

#-----------------------------------------------------------------------------
# WARNING: Don't touch anything below this line unless you exactly know it !!!
#-----------------------------------------------------------------------------


#
# Dependancies...
#
$(EXE): $(OBJ)
	$(FC) -o $@ $(FCFLAGS) $(LDFLAGS) $(OBJ) $(LDLIBS) -I$(NETCDF_INC) -L$(NETCDF_LIB) -lnetcdf -lnetcdff


modal_aero_data.o: \
	radconstants.o

gen_modal_optics_file.o: \
	radconstants.o          \
	shr_kind_mod.o          \
	module_optics_support.o \
	modal_aero_data.o       \
	miesubs.o               \
	physconst.o     

radconstants.o: \
	module_optics_support.o

shr_const_mod.o: \
	shr_kind_mod.o

physconst.o: \
	shr_kind_mod.o \
	shr_const_mod.o

clean:
	/bin/rm  *.o  *.mod *.exe *.f90 fort.*


cleanmod:
	/bin/rm  *.mod  module*.o  *.exe


