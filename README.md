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


# Fall 2021 Phys 659 Work
## WARPX




## Matrix Analysis of Single-Shot Wakefield Measurement System

This analysis is based off of the article Single-shot wakefield measurement system in the Physical Review Accelerators and Beams 21, 062801 (2018). The system outputs the beam to go through a transverse deflecting cavity (TDC), two quads (not taken into account here), and a spectrometer. 

The TDC has a transport matrix in <img src="https://latex.codecogs.com/svg.image?(x,&space;x',&space;y,&space;y',&space;z,&space;\delta)" title="(x, x', y, y', z, \delta)" /> phase space of:

<img src="https://latex.codecogs.com/svg.image?\begin{pmatrix}1&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;0\\0&space;&&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;0\\0&space;&&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;&&space;0\\0&space;&&space;0&space;&&space;0&space;&&space;1&space;&&space;\kappa&space;&&space;0\\0&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;&&space;0\\0&space;&&space;0&space;&&space;\kappa&space;&&space;0&space;&&space;0&space;&&space;1\end{pmatrix}" title="\begin{pmatrix}1 & 0 & 0 & 0 & 0 & 0\\0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & 1 & 0 & 0 & 0\\0 & 0 & 0 & 1 & \kappa & 0\\0 & 0 & 0 & 0 & 1 & 0\\0 & 0 & \kappa & 0 & 0 & 1\end{pmatrix}" />

The drift transport matrix is: 

<img src="https://latex.codecogs.com/svg.image?\begin{pmatrix}1&space;&&space;L&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\0&space;&&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\0&space;&&space;0&space;&&space;1&space;&&space;L&space;&&space;0&space;&&space;0\\0&space;&&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;&&space;0\\0&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;&&space;0\\0&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;1\end{pmatrix}" title="\begin{pmatrix}1 & L & 0 & 0 & 0 & 0 \\0 & 1 & 0 & 0 & 0 & 0 \\0 & 0 & 1 & L & 0 & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & 0 & 0 & 1 & 0\\0 & 0 & 0 & 0 & 0 & 1\end{pmatrix}" />

And, finally, the spectrometer acts like a dipole and has a transport matrix of: 

<img src="https://latex.codecogs.com/svg.image?\begin{pmatrix}\cos(\theta)&space;&&space;\rho&space;\sin(\theta)&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;\rho&space;(1-\cos(\theta))\\-\sin(\theta)/\rho&space;&&space;\cos(\theta)&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;\sin(\theta)\\0&space;&&space;0&space;&&space;1&space;&&space;\rho\theta&space;&&space;0&space;&&space;0\\0&space;&&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;&&space;0\\-\sin(\theta)&space;&&space;-\rho&space;(1-\cos(\theta))&space;&&space;0&space;&&space;0&space;&&space;1&space;&&space;\rho(\sin(\theta)&space;-&space;\theta)\\0&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;&&space;1\end{pmatrix}" title="\begin{pmatrix}\cos(\theta) & \rho \sin(\theta) & 0 & 0 & 0 & \rho (1-\cos(\theta))\\-\sin(\theta)/\rho & \cos(\theta) & 0 & 0 & 0 & \sin(\theta)\\0 & 0 & 1 & \rho\theta & 0 & 0\\0 & 0 & 0 & 1 & 0 & 0\\-\sin(\theta) & -\rho (1-\cos(\theta)) & 0 & 0 & 1 & \rho(\sin(\theta) - \theta)\\0 & 0 & 0 & 0 & 0 & 1\end{pmatrix}" />

The code for this is a simple jupyter notebook, but does require a specific library of sobol. For more information, go to <a href="https://people.sc.fsu.edu/~jburkardt/py_src/sobol/sobol.html">Sobol</a>. 

The single jupyter notebook is all that is needed and will be commented out shortly.


