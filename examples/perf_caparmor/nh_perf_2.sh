#!/bin/bash

# Number of procs in power of 2
nbpi=3
nbpf=4
#nbpf=4

# Horizontal domain resolution in power of 2
resi=11
resf=11
#resf=11

# Test name
tn='mainmodel'
sd='/home2/caparmor/grima/NH_Multigrid/IntelMPI/bin'

for (( n=nbpi; n<=nbpf; n++ ))
do

   echo""
   nbp=$(( 2 ** $n ))

  for (( r=resi; r<=resf; r++ ))
  do
 
    res=$(( 2 ** $r ))

    echo "${nbp}x${nbp}_${res}:"
    echo "-------"

    dir="${nbp}x${nbp}_${res}"
    echo "- Making directory: $dir"
    \mkdir $dir
    cd $dir

    nml="namelist_${nbp}x${nbp}_${res}"
    echo "- creating namelist file: $nml"
    sed 's/RES/'$res'/g' < ../namelist_ref > $nml
    sed -i 's/NBP/'$nbp'/g' $nml
    ln -fs $nml 'namelist'

    echo "- Linking with executable test: $tn"
    ln -fs ${sd}/${tn} ${tn}

    cd ..
  done 

  nn=$(( $nbp  * $nbp ))
  init=$(( 2 ** $resi ))
  final=$(( 2 ** $resf ))
  pbsr="run_${nn}p.pbs"
  pbs="run_${nn}p_${tn}.pbs"
  echo "- Creating PBS file for tests with $n procs: ${pbs}"
  \rm ${pbs}
  sed    's/INIT/'${init}'/g' < ${pbsr} > ${pbs}
  sed -i 's/FINAL/'${final}'/g' ${pbs}
  sed -i 's/NHTEST/'$tn'/g' ${pbs}
  echo "="
  echo "=="
  echo "===> Submmit ${pbs} "
  echo "=="
  echo "="
  qsub ${pbs}

done
