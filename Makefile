
MAKEFILE_IN = $(PWD)/Makefile.in
include $(MAKEFILE_IN)

MAKEFILE_MACH  = $(PWD)/Makefile.mach
include $(MAKEFILE_MACH)

TEMPLATES = bexp bx3d earth fieldloop flock flock_grmhd isentropic jupiter kepler kh mri2 rotor shear shocktube spinring spread vortex sorathia_grmhd blast_grmhd fieldloop_grmhd bl kep_ring torus_fm isentropic_RAM entropywave acousticwave cb acousticwave_cart bondi binarybondi alfvenwave alfvenwave_cart magnetosonicwave_cart magnetosonicwave_cart3d bondi_mhd fieldloop_linear_cart ecc tde acousticwave_sph cartesian_shear cartesian_shear_sph vortexShedding

GIT_VERSION = $(shell git describe --dirty --always --tags)

OPT_DEFS = -DGIT_VERSION=\"$(GIT_VERSION)\"
OPT_DEFS += -DNUM_C=$(NUM_C)
OPT_DEFS += -DNUM_N=$(NUM_N)
OPT_DEFS += -DCT_MODE=$(CT_MODE)
OPT_DEFS += -DTYPE=$(TYPE)
OPT_DEFS += -DSUBTYPE=$(SUBTYPE)

FLAGS = -O3 -Wall -g $(OPT_DEFS) $(DIR_DEFS)

INC = -I$(H55)/include
LIB = -L$(H55)/lib -lhdf5 -lm

OBJ = main.o readpar.o faces.o grid.o profiler.o geometry.o geometry_functions.o hydro.o timestep.o domain_aos.o substep_aos.o plm_aos.o domain_soa1.o substep_soa1.o plm_soa1.o

default: disco

.PHONY: $(TEMPLATES)

$(TEMPLATES):
	cp Templates/$@.par in.par
	cp Templates/$@.in Makefile.in
	make clean
	make

%.o: %.c header.h
	$(CC) $(FLAGS) $(LOCAL_CFLAGS) $(INC) -c $< -o $@

disco: $(OBJ) header.h
	$(CC) $(FLAGS) -o disco $(OBJ) $(LOCAL_LDFLAGS) $(LIB)

clean:
	rm -f *.o disco
