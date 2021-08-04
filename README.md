# Summer 2021 Phys 659 Work

These codes are to simulate a Structure-Based Wakefield Acceleration. 

## Tremaine-Roz C++ Code

A C++ based simulation code that has two possibilities for runs: an on-axis representation of what's going on (including E and B fields) and a scan over x and y at some target z. 

### General Comments
#### Files to keep together
In order to run the main code ('run_tr_3d_scan.cpp' or 'run_tr_3d.cpp'), two other files are required to be in the same directory, which are the 'tremaine_roz.cpp' and 'tremaine_roz.h'. 

#### How to run this code
Once the three files are together in a directory, the Windows commands are:
'g++ mainfile.cpp -lm -o wake' - This looks through the three files and creates an executable file
'.\wake.exe' - This executable file actually runs the code and, for the scan, will have output as the scan processes the data for each (x,y) point. 

### On-Axis run
'run_tr_3d.cpp' is for the on-axis run and '3dscan' is for the scan over x and y. The output 



### Scan over X and Y
The output is 
