#!/bin/bash
i=0;
max=60;
n=1000000;
while [ "$i" -le "$max" ]; do
  mkdir "s17E_$i";
  cd "s17E_$i";
  skip=$(python -c "print 0+$n*$i");
  echo $'#!/bin/sh' > script.sh
  echo $'#BSUB -o test.log' >> script.sh
  echo $'eval `scram runtime -sh`' >> script.sh
  echo "pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/charmonium2017Lists/Charmonium_Run2017E_DCAP.list hist$i.root -v outputFile ntu$i.root -v histoMode RECREATE -v use_gen f -v process BuJPsiK -n $n -s $skip" >> script.sh
  echo "" >> script.sh
  bsub < script.sh;
  cd ..;
  i=`expr "$i" + 1`;
done
