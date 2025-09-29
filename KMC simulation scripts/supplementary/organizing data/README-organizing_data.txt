To convert the simulation output data, which stores vertex coordinates as lammps .dat files and the faces as .vtk files, we perform the following. 
1. Run grabCoordinates.sh on the high performace computing cluster (HPCC), to convert the vertex coordinate .dat lammps files into CSVs. 
2. Run grabFaces.sh on the HPCC to convert the face .vtk files into CSVs.
3. Transfer the coordinate and face data to a local machine, (e.g., using WINSCP).
4. Run convertSpaces2Commas.py to enfore the CSV file format, for both the coordinate and face data. This completes the process of converting the output data into the following standard format:
 - Coordinates: Nx3 matrix of x,y,z coordinates, where each row is a vertex, i.e., [X Y Z]
 - Faces: Mx4 matrix, where each row is a face, columns 1-3 are vertex indices of from the corresponding coordinates matrix, and column 4 is the particle type of that triangle, i.e., [v1 v2 v3 type]