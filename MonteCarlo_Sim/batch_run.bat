@echo off
g++ -O2 MC_sim.cpp -o MC_sim.exe
start /b MC_sim.exe 100000 5.000000000000001 0
start /b MC_sim.exe 100000 17.69617959256636 1
start /b MC_sim.exe 100000 62.63095443447242 2
start /b MC_sim.exe 10000 221.66572354525292 3
start /b MC_sim.exe 10000 784.5272906745923 4
start /b MC_sim.exe 10000 2776.6271662094196 5
start /b MC_sim.exe 10000 9827.1385989681 6
start /b MC_sim.exe 10000 34780.5619056761 7
start /b MC_sim.exe 1000 123096.61396264327 8
start /b MC_sim.exe 1000 435667.95758394944 9
start /b MC_sim.exe 100 1541931.6840264306 10
start /b MC_sim.exe 100 5457260.000000001 11
