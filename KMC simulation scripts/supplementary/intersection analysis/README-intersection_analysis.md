After running the simulations, we analyze face-face intersections as follows.
1. Run convertSnapshots.sh on the high performance computing cluster (HPCC) to generate the vertex lammps .dat files and face .vtk files.
2. Transfer the resultant files to a local machine (e.g., using WINSCP).
3. Run convertCoordinates2CSV.py to convert the vertex coordinate data into CSVs.
4. Run convertVtk2Dat.sh to convert the face .vtk files into .dat files. 
5. Run ConvertFaces2CSV.py to convert the face .dat files into CSV files. 
6. Using the coordinate and face files, run countIntersectionsOverTime.py, which will save the analysis data to a JSON file. 
7. Run removeOvergrownSims.py to remove the overgrown simulation outcomes by editing the JSON file.

7. Run plotCoutsOverTime.py using the JSON data file to plot the number of intersecting faces over time. 
