#!/bin/bash

cd src;
make clean;
echo " clean "
make ham;
cd ..

for eps in  0 1 2 5
do

for delta in 0 1 2 5
do

dirout='eps_'$eps'_delta_'$delta

mkdir -p $dirout'/data'; cd $dirout

cat > INPUT_eps << EOF
${eps}
EOF

cat > INPUT_delta << EOF
${delta}
EOF

nohup ../src/ham++ > err&

echo " -- Calculating for eps = " ${eps}
echo " -- Calculating for delta = " ${delta}

cd ..

done
done
