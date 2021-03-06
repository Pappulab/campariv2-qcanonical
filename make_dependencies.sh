#!/bin/bash
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


SRC_DIR=${PWD}

LIB_DIR=${PWD}"/../lib"

MODULES="accept aminos atoms clusters commline contacts cutoffs diffrac dipolavg distrest dssps ems energies ewalds forces fos fyoc grandensembles grids inter interfaces ionize iounit keys martin_own math mcgrid mcsums mini molecule movesets mpistuff paircorr params pdb polyavg polypep sequen shakeetal system tabpot threads torsn ujglobals units wl zmatrix"

PRESOURCES="accsim allocate assignprm backup boundary cartld cartmd chainsaw clustering clustering_utils conrot constraint_solvers dssp emicroscopy energy energy_wrap ensemble ewald flow fmcscgrid fmsmcmpi force force_wrap fyczmat getkey graph_algorithms holes initial inner_loops inner_loops_en inner_loops_imp intld intmd makeio makepept math_utils mcmove mcstat minimize nucconrot parsefiles parsekey particlefluctuation polar polymer proteus prtpdb readfyc readgrid readpdb readprm restart rigidmoves sanity_checks sav sav_auxil setconf sidechain string_utils structure summary titrate topology torconrot torsion ujconrot ujsugar_pucker unbond utilities wanglandau"

SOURCES=""
for i in $PRESOURCES;
do
  SOURCES=${SOURCES}"${i}.f90 "
done

MODOS=""
for i in $MODULES;
do
  MODOS=${MODOS}"mod_${i}.f90 "
done

rm ${SRC_DIR}/DEPENDENCIES
for i in $MODULES;
do
  echo \${LIB_DIR}/\${ARCH}/${i}.o \${LIB_DIR}/\${ARCH}/mpi/${i}.o \${LIB_DIR}/\${ARCH}/${i}.mod \${LIB_DIR}/\${ARCH}/mpi/${i}.mod: \${SRC_DIR}/mod_${i}.f90 >> ${SRC_DIR}/DEPENDENCIES
  DEPS=`grep -H -m 1 "  use ${i}" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/" $1}' | sed s/.f90:/.o/`
  DEPS2=`grep -H -m 1 "  use ${i}" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/mpi/" $1}' | sed s/.f90:/.o/`
  DEPS3=`grep -H -m 1 "  use ${i}" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/" $1}' | sed s/.f90:/.o/`
  DEPS4=`grep -H -m 1 "  use ${i}" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/mpi/" $1}' | sed s/.f90:/.o/`
  echo ${DEPS3} ${DEPS}: \${LIB_DIR}/\${ARCH}/${i}.o \${LIB_DIR}/\${ARCH}/${i}.mod >> ${SRC_DIR}/DEPENDENCIES
  echo ${DEPS4} ${DEPS2}: \${LIB_DIR}/\${ARCH}/mpi/${i}.o \${LIB_DIR}/\${ARCH}/mpi/${i}.mod >> ${SRC_DIR}/DEPENDENCIES
done

for i in $SOURCES;
do
  echo \${LIB_DIR}/\${ARCH}/${i%.f90}.o \${LIB_DIR}/\${ARCH}/mpi/${i%.f90}.o: \${SRC_DIR}/${i} >> ${SRC_DIR}/DEPENDENCIES
done

MDEPS=`grep -H -m 1 "#include \"macros.i\"" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/" $1}' | sed s/.f90:#include/.o/`
MDEPS2=`grep -H -m 1 "#include \"macros.i\"" ${SOURCES} | awk '{print "${LIB_DIR}/${ARCH}/mpi/" $1}' | sed s/.f90:#include/.o/`
echo ${MDEPS} ${MDEPS2}: \${SRC_DIR}/macros.i >> ${SRC_DIR}/DEPENDENCIES

MDEPS=`grep -H -m 1 "#include \"macros.i\"" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/" $1}' | sed s/.i:#include/.o/`
MDEPS2=`grep -H -m 1 "#include \"macros.i\"" ${MODOS} | awk '{print "${LIB_DIR}/${ARCH}/mpi/" $1}' | sed s/.i:#include/.o/`
echo ${MDEPS} ${MDEPS2}: \${SRC_DIR}/macros.i >> ${SRC_DIR}/DEPENDENCIES

