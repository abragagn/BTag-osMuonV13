#!/bin/bash
i=0;
max=24;
n=1000000;
while [ "$i" -le "$max" ]; do
  mkdir "s17_$i";
  cd "s17_$i";
  skip=$(python -c "print 0+$n*$i");
  echo $'#!/bin/sh' > script.sh
  echo $'#BSUB -o test.log' >> script.sh
  echo $'eval `scram runtime -sh`' >> script.sh
  echo "pdTreeAnalyze /lustre/cmswork/abragagn/ntuList/MC2017Lists/BuToJpsiK_2017_DCAP.list hist$i.root -v outputFile ntu$i.root -v histoMode RECREATE -v use_gen t -v process BuJPsiK -n $n -s $skip" >> script.sh
  echo "" >> script.sh
  bsub < script.sh;
  cd ..;
  i=`expr "$i" + 1`;
done
