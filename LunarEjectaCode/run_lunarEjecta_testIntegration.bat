@echo off
del a.exe lunarEjecta_TestIntegration.o lunarEjecta_FractalIntegration.o lunarEjecta_GeneralExpressions.o lunarEjecta_AdaptiveMesh.o
g++ -std=c++11 -O2 -c lunarEjecta_TestIntegration.cpp lunarEjecta_FractalIntegration.cpp lunarEjecta_GeneralExpressions.cpp lunarEjecta_AdaptiveMesh.cpp
g++ lunarEjecta_TestIntegration.o lunarEjecta_FractalIntegration.o lunarEjecta_GeneralExpressions.o lunarEjecta_AdaptiveMesh.o
a.exe