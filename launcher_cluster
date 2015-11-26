#!/bin/bash

rm Beta_* -rf

cd src;
make wt;
cd ..

echo "Compilation done "

rm temperatures
cat >> temperatures << EOF
EOF

for beta in 10 13 16 20 25 32 40 50 63 80 100
do

printf ${beta} >> temperatures
printf " " >> temperatures

dirout='Beta_'$beta

mkdir -p $dirout'/data'; cd $dirout

cat > input_file << EOF
Nombre de repetitions       : 100000
Nombre de temps de sortie   : 1
Magnitude deplacement       : 0.1
Method                      : 0
1/(K x Temperature)         : ${beta}
Nombre de dimensions        : 2
Frequence acquisition       : 100
Borne zone de biais         : 1.2
Pas d'espace                : 0.1
Exposant WL                 : 0.8
Constante WL                : 1.
Parametre Well-Tempered     : 1.0
EOF

nohup ../src/WT++ > err&

echo " -- Calculating for beta = " ${beta}

cd ..

done

printf " \n" >> temperatures 
