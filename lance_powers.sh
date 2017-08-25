#!/bin/bash

cd src;
make clean;
echo " clean "
make ham;
cd ..

for eps in  0.5 0.8 1.0 1.2 1.5 1.8 2.0 2.2 2.5 2.8 3.0 4.0
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

cd ..


done
