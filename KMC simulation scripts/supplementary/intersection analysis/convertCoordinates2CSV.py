import os
import numpy as np
import sys
from pathlib import Path

def convertCoordinates2CSV(inputDirectory, outputDirectory):
    coordinatesStartLine = 14
    indxEnd = 2

    dirList = os.listdir(inputDirectory)
    for folderName in dirList:
        print(f"Processing {folderName}")
        newFolderPath = os.path.join(outputDirectory, folderName)
        if not os.path.exists(newFolderPath):
            os.mkdir(newFolderPath)
        fileList = os.listdir(os.path.join(inputDirectory, folderName, 'conversions'))
        for fileName in fileList:
            if "lammps" in fileName:
                coordinates = []
                with open(os.path.join(inputDirectory, folderName, 'conversions', fileName), "r") as file:
                    for (i, line) in enumerate(file):
                        if i >= coordinatesStartLine:
                            if  "Bonds" in line:
                                break
                            if line!="\n":
                                lineVals = [float(x) for x in line.split(" ") if x!="\n"]
                                coordinates.append(lineVals[indxEnd:])
                coordinates = np.array(coordinates)
                newFileName = fileName.replace(".dat", "_coordinates.CSV")
                np.savetxt(os.path.join(newFolderPath, newFileName), coordinates, delimiter=",")

def main():
    # Usage: python convertVtk2CSV.py "path/to/input_directory" "path/to/output_directory"
    userArgs = sys.argv
    inputDirectory = Path(userArgs[1])
    outputDirectory = Path(userArgs[2])

    convertCoordinates2CSV(inputDirectory,outputDirectory)

if __name__ == '__main__':
    main()