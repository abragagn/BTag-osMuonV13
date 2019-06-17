#!/bin/bash
i=0;
max=101;
echo $'#!/bin/sh' > park.sh
echo -n "hadd ntuBsDG0MC2018.root" >> park.sh
while [ "$i" -le "$max" ]; do
  echo -n " s$i/ntu$i" >> park.sh
  echo -n ".root" >> park.sh
  i=`expr "$i" + 1`;
done
echo " " >> park.sh
bash park.sh;
rm park.sh;
mv ntuBsDG0MC2018.root ../
