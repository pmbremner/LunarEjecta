@echo off
g++ -O2 MC_sim.cpp -o MC_sim.exe
MC_sim.exe 10000 1000. 0
python plotTraj.py