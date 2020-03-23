@echo off
del a.exe lunarEjecta_Main.o lunarEjecta_MeteoroidFlux.o 
g++ -std=c++11 -O2 -c lunarEjecta_Main.cpp lunarEjecta_MeteoroidFlux.cpp
g++ lunarEjecta_Main.o lunarEjecta_MeteoroidFlux.o
a.exe