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

### Running the code
This code is run on the NICADD cluster as there is no `WARP` option for Windows and because the HDF5 files are rather memory intensive. It's python-based, so a connection between the cluster and the personal Windows computer is created to run Jupyterlab through Chrome. A Diag directory is created and the HDF5 files dumped for certain timesteps (defined in the file). 

### HDF5 Files and Plot Creation
The HDF5 file has a nested format for the data and generally needs investigating if unfamiliar with the output of a specific run. Below is the current layoutfor the HDF5 files created from the `testingawacode-WORKS.ipynb`. 

```
-Data
 -Timestep (eg "500")
  -Fields
   -B
    -x (data)
    -y (data)
    -z (data)
   -E
    -x (data)
    -y (data)
    -z (data)
   -J
    -x (data)
    -y (data)
    -z (data)
   -Rho (data)
  -Particles
   -Helium0+
    -Momentum 
     -x (data)
     -y (data)
     -z (data)
    -positionOffset
     -x (data)
     -y (data)
     -z (data)
   ...
```


Plotting is done using `WARPplotting.ipynb`, where the HDF5 file is read and data taken from it for a specific timestep, described at the beginning of the file.




