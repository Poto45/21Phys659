# Summer 2021 Phys 659 Work

These codes are used to simulate a Structure-Based Wakefield Acceleration. 

## Tremaine-Roz C++ Code

A C++ based simulation code that has two possibilities for runs: an on-axis representation of what's going on (including E and B fields) and a scan over x and y at some target z. 

### General Comments
#### Files to keep together
In order to run the main code (`run_tr_3d_scan.cpp` or `run_tr_3d.cpp`), two other files are required to be in the same directory, which are the `tremaine_roz.cpp` and `tremaine_roz.h`. 

#### How to run this code
Once the three files are together in a directory, the Windows commands are:

`g++ mainfile.cpp -lm -o wake` - This looks through the three files and creates an executable file

`.\wake.exe` - This executable file actually runs the code and, for the scan, will have output as the scan processes the data for each (x,y) point. 

### On-Axis run
`run_tr_3d.cpp` is for the on-axis run. The output is given in two files: `wake_ExEyEzBxByBz.dat` and `wake_z_wz.dat`. 



### Scan over X and Y
The output electric and magnetic field data is given in the file `wake_2Dsurf.dat` with extra data given in `wake_z_wz.dat`. Other output can also be chosen, but see the code in here and documentation to add other parameters. The `wake_2dsurf.dat` file is what has the x, y, z, Ex, Ey, Bx, and By data and are delimited using a space and can be extracted using python to plot (`ScannedPlots.ipynb`). The layout is seen below:

```
  x      y    z   Ex  Ey  Bx  By
-----  ----- --- --- --- --- ---
 -5.5   1.5   0   0   0   0   0
...
 5.5    1.5   0   0   0   0   0
```


### Plotting
The wakefields are plotted in the file `ScannedPlots.ipynb`. In order to get the wakefields, the equation ![equation](https://latex.codecogs.com/gif.latex?F%28x%2Cy%29%20%3D%20%28E_x%20-%20cB_y%2C%20E_y%20&plus;%20cB_x%29%5ET) is implemented. 





## WARP Code
This python-based code is powerful in dumping large amounts of data into HDF5 files and continuing its calculations through many different timesteps. 







