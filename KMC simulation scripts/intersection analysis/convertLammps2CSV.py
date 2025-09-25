import os
import numpy as np
import argparse

def main():
    coordinatesStartLine = 14
    indxEnd = 2

    inputDirectory = r"C:\Users\mason\OneDrive\Desktop\code\openMesh\test_20250924_toroid\highResSamples\on_target_R_exc0_15"
    outputDirectory = r"C:\Users\mason\OneDrive\Desktop\code\openMesh\test_20250924_toroid\highResSamples\on_target_R_exc0_15_coordinates"
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


if __name__ == '__main__':
    main()