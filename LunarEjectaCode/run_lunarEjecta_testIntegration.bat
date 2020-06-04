@echo off
del a.exe lunarEjecta_TestIntegration.o lunarEjecta_FractalIntegration.o lunarEjecta_GeneralExpressions.o
g++ -std=c++11 -O2 -c lunarEjecta_TestIntegration.cpp lunarEjecta_FractalIntegration.cpp lunarEjecta_GeneralExpressions.cpp
g++ lunarEjecta_TestIntegration.o lunarEjecta_FractalIntegration.o lunarEjecta_GeneralExpressions.o
a.exe