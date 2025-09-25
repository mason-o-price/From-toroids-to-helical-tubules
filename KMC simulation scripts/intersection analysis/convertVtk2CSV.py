import os
import numpy as np
import argparse

def main():
    indxEnd = 1

    inputDirectory = r"C:\Users\mason\OneDrive\Desktop\code\openMesh\test_20250924_toroid\highResSamples\on_target_R_exc0_15"
    outputDirectory = r"C:\Users\mason\OneDrive\Desktop\code\openMesh\test_20250924_toroid\highResSamples\on_target_R_exc0_15_faces"
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

if __name__ == '__main__':
    main()