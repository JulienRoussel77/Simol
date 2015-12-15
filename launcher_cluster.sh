#!/bin/bash

for nbParticles in 64 128 256 512 1024 2048 4096 8192
do

mkdir -p simol/output/BoundaryLangevin/TriChain/Quadratic/N${nbParticles}

rm -f cluster_input_trichain_N${nbParticles}.yaml
cat > cluster_input_trichain_N${nbParticles}.yaml << EOF
Geometry:
  Dimension: 1
  Length: 1

Mesh:
  Time:
    Step: .01
    FinalTime: 100

Physics:
  System:
    Name: TriChain
    Number: ${nbParticles}
    Position: [0, 0]
  Potential:
    Name: Quadratic
    Alpha: 2
  Model:
    Name: BoundaryLangevin
    Gamma: 1
    TemperatureLeft: 1.1
    TemperatureRight: .9
    Tau: 0

Output:
  Foldername: N${nbParticles}
  FinalFile: 
  DecorrelationTime: 2
  Period: 10

ControlVariate:
  Name: None
EOF

cd simol/build
nohup src/molecular_dynamics -i "../../cluster_input_trichain_N${nbParticles}.yaml" > ../output/BoundaryLangevin/TriChain/Quadratic/cluster_N${nbParticles}.out&
cd ../..

echo " -- Calculating for N = " ${nbParticles}

done
