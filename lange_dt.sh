#!/bin/bash

cd src;
make clean;
echo " clean "
make ham;
cd ..

for eps in  0 2 6
do

for dt in 0.005 0.01 0.012 0.015
do

dirout='dt_'$dt'eps_'$eps

mkdir -p $dirout'/data'; cd $dirout

cat > INPUT_dt << EOF
${dt}
EOF

cat > INPUT_eps << EOF
${eps}
EOF


nohup ../src/ham++ > err&

echo " -- Calculating for dt = " ${dt}
echo " -- Calculating for eps = " ${eps}

cd ..


done
done
