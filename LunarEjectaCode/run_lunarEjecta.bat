@echo off
del a.exe lunarEjecta_Main.o lunarEjecta_GeneralExpressions.o lunarEjecta_MeteoroidFlux.o lunarEjecta_SecondaryFluxData.o lunarEjecta_Assembly.o
g++ -std=c++11 -O2 -c lunarEjecta_Main.cpp lunarEjecta_GeneralExpressions.cpp lunarEjecta_MeteoroidFlux.cpp lunarEjecta_SecondaryFluxData.cpp lunarEjecta_Assembly.cpp
g++ lunarEjecta_Main.o lunarEjecta_GeneralExpressions.o lunarEjecta_MeteoroidFlux.o lunarEjecta_SecondaryFluxData.o lunarEjecta_Assembly.o
a.exe