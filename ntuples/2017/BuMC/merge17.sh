#!/bin/bash
i=0;
max=24;
echo $'#!/bin/sh' > park.sh
echo -n "hadd ntuBuMC2017.root" >> park.sh
while [ "$i" -le "$max" ]; do
  echo -n " s17_$i/ntu$i" >> park.sh
  echo -n ".root" >> park.sh
  i=`expr "$i" + 1`;
done
echo " " >> park.sh
bash park.sh;
rm park.sh;
mv ntuBuMC2017.root ../