#!/bin/sh
constraint=$1
chunk=$2

make -j
./ParasoR --pre --constraint $constraint --stem -r > ../doc/pre.txt
echo "" > ../doc/stem.txt
for i in `seq 0 ${chunk}`
do
    ./ParasoR --constraint $constraint -i $i -k $chunk -r >> ../doc/stem.txt
done
./ParasoR --constraint $constraint --stem -r >> ../doc/stem.txt

echo "" > ../doc/stem_mem.txt
for i in `seq 0 ${chunk}`
do
    ./ParasoR --save_memory --constraint $constraint -i $i -k $chunk -r >> ../doc/stem_mem.txt
done
for i in `seq 0 ${chunk}`
do
    ./ParasoR --save_memory --constraint $constraint -i $i -k $chunk -r >> ../doc/stem_mem.txt
done
./ParasoR --constraint $constraint --stem  -r >> ../doc/stem_mem.txt
