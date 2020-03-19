@echo off
del a.exe lunarEjecta_Main.o lunarEjecta_MeteoroidFlux.o vector2d.o
g++ -std=c++11 -O2 -c lunarEjecta_Main.cpp lunarEjecta_MeteoroidFlux.cpp vector2d.cpp
g++ lunarEjecta_Main.o lunarEjecta_MeteoroidFlux.o vector2d.o
a.exe