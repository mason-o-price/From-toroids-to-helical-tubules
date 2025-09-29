import os
import numpy as np
import sys
from pathlib import Path

def convertFaces2CSV(inputDirectory, outputDirectory):
    indxEnd = 1
    dirList = os.listdir(inputDirectory)
    for folderName in dirList:
        print(f"Processing {folderName}")
        newFolderPath = os.path.join(outputDirectory, folderName)
        if not os.path.exists(newFolderPath):
            os.mkdir(newFolderPath)
        fileList = os.listdir(os.path.join(inputDirectory, folderName, 'conversions'))
        for fileName in fileList:
            if "vtk" in fileName and ".dat" in fileName:
                faces = []
                with open(os.path.join(inputDirectory, folderName, 'conversions', fileName), "r") as file:
                    faceLineStart = None
                    faceLineEnd = None
                    # Find the start of the face data
                    for (lineIdx, line) in enumerate(file):
                        if "CELLS" in line:
                            faceLineStart = lineIdx
                        if "CELL_TYPES" in line:
                            faceLineEnd = lineIdx
                            break
                    # Read the face data
                with open(os.path.join(inputDirectory, folderName, 'conversions', fileName), "r") as file:
                    for (lineIdx, line) in enumerate(file):
                        if faceLineStart is None or faceLineEnd is None:
                            Warning(f"Failed to find face line start for file {fileName}")
                        if lineIdx > faceLineStart and lineIdx < faceLineEnd:
                            lineVals = [int(x) for x in line.split(" ") if x!="\n"]
                            faces.append(lineVals[indxEnd:])
                faces = np.array(faces)
                # print(faces)
                newFileName = fileName.replace(".dat", "_faces.CSV")
                np.savetxt(os.path.join(newFolderPath, newFileName), faces, delimiter=",")

def main():
    # Usage: python convertVtk2CSV.py "path/to/input_directory" "path/to/output_directory"
    userArgs = sys.argv
    inputDirectory = Path(userArgs[1])
    outputDirectory = Path(userArgs[2])

    convertFaces2CSV(inputDirectory,outputDirectory)
    

if __name__ == '__main__':
    main()