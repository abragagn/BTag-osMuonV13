#!/bin/bash
i=0;
max=212;
echo $'#!/bin/sh' > park.sh
echo -n "hadd ntuBuData2018D.root" >> park.sh
while [ "$i" -le "$max" ]; do
  echo -n " s18D_$i/ntu$i" >> park.sh
  echo -n ".root" >> park.sh
  i=`expr "$i" + 1`;
done
echo " " >> park.sh
bash park.sh;
rm park.sh;
mv ntuBuData2018D.root ../
