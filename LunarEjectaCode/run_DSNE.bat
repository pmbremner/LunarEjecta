@echo off
.\a.exe .\param_files\parameters_equator_DSNE.txt run_equator_DSNE_run0.txt > test_equator_DSNE_run0.txt
.\a.exe .\param_files\parameters_lat45_DSNE.txt run_lat45_DSNE_run0.txt > test_lat45_DSNE_run0.txt
.\a.exe .\param_files\parameters_pole_DSNE.txt run_equator_DSNE_run0.txt > test_pole_DSNE_run0.txt