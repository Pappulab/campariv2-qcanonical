#--------------------------------------------------------------------------#
# LICENSE INFO:                                                            #
#--------------------------------------------------------------------------#
#    This file is part of CAMPARI.                                         #
#                                                                          #
#    Version 2.0                                                           #
#                                                                          #
#    Copyright (C) 2014, The CAMPARI development team (current and former  #
#                        contributors)                                     #
#                        Andreas Vitalis, Adam Steffen, Rohit Pappu, Hoang #
#                        Tran, Albert Mao, Xiaoling Wang, Jose Pulido,     #
#                        Nicholas Lyle, Nicolas Bloechliger                #
#                                                                          #
#    Website: http://sourceforge.net/projects/campari/                     #
#                                                                          #
#    CAMPARI is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by  #
#    the Free Software Foundation, either version 3 of the License, or     #
#    (at your option) any later version.                                   #
#                                                                          #
#    CAMPARI is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
#    GNU General Public License for more details.                          #
#                                                                          #
#    You should have received a copy of the GNU General Public License     #
#    along with CAMPARI.  If not, see <http://www.gnu.org/licenses/>.      #
#--------------------------------------------------------------------------#
# AUTHORSHIP INFO:                                                         #
#--------------------------------------------------------------------------#
#                                                                          #
# MAIN AUTHOR:   Andreas Vitalis                                           #
# CONTRIBUTIONS: Adam Steffen, Albert Mao                                  #
#                                                                          #
#--------------------------------------------------------------------------#

CAMPARI_HOME=/packages/campari
                                                                                                                                  
# x86_64 locale (Default)
ARCH=x86_64
BIN_DIR=${CAMPARI_HOME}/bin
LIB_DIR=${CAMPARI_HOME}/lib
SRC_DIR=${CAMPARI_HOME}/source
FF=gfortran
LLFLAGS=-O3
EXTRA_LIBS=
MPIFF=mpif90
MPILFLAGS=-O2
MPIEXTRA_LIBS=
MV=/bin/mv
RM=/bin/rm
INTELDEFAULTS=-O3 -r8 -qoverride_limits -stand f03 -fpp
INTELDEBUG=-check bounds -check uninit -check arg_temp_created -check pointers -traceback -g -debug extended -debug-parameters all
SUNDEFAULTS=-r8const -fpp -xO4 -ftrap=common,no%overflow
SUNDEBUG=-xcheck=init_local -g
GNUDEFAULTS=-cpp -O3 -fall-intrinsics -ffpe-trap=invalid,zero -std=f2003 -fdefault-real-8 -fdefault-double-8
GNUDEBUG=-fbacktrace -fbounds-check -fcheck-array-temporaries -g
COMPDEFAULTS=${GNUDEFAULTS} # ${GNUDEBUG}

# Please define any local settings in the file Makefile.local and avoid certain strings in pathnames, as they
# are run through patsubst (see below)
include Makefile.local

# Include the script-generated module dependencies: note this uses LIB_DIR, SRC_DIR and ARCH
include DEPENDENCIES

# to add a module, add the module here, add the module in make_dependencies.sh, re-run make_dependencies, and compile
PREPREMODS=accept aminos atoms clusters commline contacts cutoffs diffrac dipolavg distrest dssps ems energies ewalds forces fos fyoc grandensembles grids inter interfaces ionize iounit keys martin_own math mcgrid mcsums mini molecule movesets mpistuff paircorr params pdb polyavg polypep sequen shakeetal system tabpot threads torsn ujglobals units wl zmatrix

PREMODS=$(addsuffix .f90,$(PREPREMODS))
MODULES=$(patsubst %.f90,$(SRC_DIR)/mod_%.f90,$(PREMODS)) 
MODOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/%.o,$(PREMODS))
MPIMODOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/mpi/%.o,$(PREMODS))

print:
	${MODOBJS}

# to add a source, add the source here, add the source in make_dependencies.sh, re-run make_dependencies, and compile
PREPRESRCS=accsim allocate assignprm backup boundary cartld cartmd chainsaw clustering clustering_utils conrot constraint_solvers dssp emicroscopy energy energy_wrap ensemble ewald flow fmcscgrid fmsmcmpi force force_wrap fyczmat getkey graph_algorithms holes initial inner_loops inner_loops_en inner_loops_imp intmd intld makeio makepept math_utils mcmove mcstat minimize nucconrot parsefiles parsekey particlefluctuation polar polymer proteus prtpdb readfyc readgrid readpdb readprm restart rigidmoves sanity_checks sav sav_auxil setconf sidechain string_utils structure summary titrate topology torconrot torsion ujconrot ujsugar_pucker unbond utilities wanglandau

PRESRCS=$(addsuffix .f90,$(PREPRESRCS))
SOURCES=$(patsubst %.f90,${SRC_DIR}/%.f90,$(PRESRCS))
OBJECTS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/%.o,$(PRESRCS))
MPIOBJS=$(patsubst %.f90,${LIB_DIR}/${ARCH}/mpi/%.o,$(PRESRCS))

# to correct for the m2c error
%.o : %.mod


# our strategy is to keep objects completely out of the SRC-directory (not even temporary) and to use a sub-directory for the MPI-version
# that means all the targets are only defined in that target directory, which is necessary in order to be able to address, maintain, and
# remove them independently.
# note that this includes the module interface definitions (.mod), which are unfortunately not compiler-independent. but this way it at least
# allows us to do preprocessing in modules without defining separate modules
# also note that dependencies are and should be exclusively handled through DEPENDENCIES
${OBJECTS} :
	${FF} ${FFLAGS} -c ${COMPDEFAULTS} -I${LIB_DIR}/${ARCH} $(patsubst ${LIB_DIR}/${ARCH}/%.o,${SRC_DIR}/%.f90,$@) -o $@

${MPIOBJS} :
	${MPIFF} ${MPIFFLAGS} -c ${COMPDEFAULTS} -DENABLE_MPI -I${LIB_DIR}/${ARCH}/mpi $(patsubst ${LIB_DIR}/${ARCH}/mpi/%.o,${SRC_DIR}/%.f90,$@) -o $@

${MODOBJS} :
	${FF} ${FFLAGS} -c ${COMPDEFAULTS} -I${LIB_DIR}/${ARCH} $(patsubst ${LIB_DIR}/${ARCH}/%.o,${SRC_DIR}/mod_%.f90,$@) -o $@
	${MV} $(patsubst ${LIB_DIR}/${ARCH}/%.o,${SRC_DIR}/%.mod,$@) $(patsubst %.o,%.mod,$@)

${MPIMODOBJS} :
	${MPIFF} ${MPIFFLAGS} -c ${COMPDEFAULTS} -I${LIB_DIR}/${ARCH}/mpi -DENABLE_MPI $(patsubst ${LIB_DIR}/${ARCH}/mpi/%.o,${SRC_DIR}/mod_%.f90,$@) -o $@
	${MV} $(patsubst ${LIB_DIR}/${ARCH}/mpi/%.o,${SRC_DIR}/%.mod,$@) $(patsubst %.o,%.mod,$@)

# library assembly is trivial (note that the name of the target doesn't match the target -> forced execution every time!)
library: $(OBJECTS)
	ar -rclvs ${LIB_DIR}/${ARCH}/lcampari.a ${MODOBJS} ${OBJECTS}

library_mpi: $(MPIOBJS)
	ar -rclvs ${LIB_DIR}/${ARCH}/mpi/lcampari_mpi.a ${MPIMODOBJS} ${MPIOBJS}

# linking hopefully as well (see above note, true here due to target location!)
campari: library
	${FF} $(LFLAGS) -o ${BIN_DIR}/${ARCH}/campari ${LIB_DIR}/${ARCH}/chainsaw.o ${LIB_DIR}/${ARCH}/lcampari.a ${EXTRA_LIBS}

campari_mpi: library_mpi
	${MPIFF} $(MPILFLAGS) -DENABLE_MPI -o ${BIN_DIR}/${ARCH}/campari_mpi ${LIB_DIR}/${ARCH}/mpi/chainsaw.o /project/fava/previous-packages/xdr/xdrf/libxdrf.a ${LIB_DIR}/${ARCH}/mpi/lcampari_mpi.a ${MPIEXTRA_LIBS} 


# some fake targets, which are really clean-up commands
objclean:
	for i in ${LIB_DIR}/${ARCH}/*.o; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi/*.o; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/*.mod; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi/*.mod; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done
	for i in ${LIB_DIR}/${ARCH}/mpi/lcampari_mpi.a ${LIB_DIR}/${ARCH}/lcampari.a; do (if [ -e $${i} ]; then ${RM} $${i}; fi;) done

clean: objclean
	for i in ${BIN_DIR}/${ARCH}/campari*; do (if [ -e $$i ]; then ${RM} $$i; fi;) done

