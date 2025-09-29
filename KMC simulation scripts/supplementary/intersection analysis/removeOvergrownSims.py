import os
import json
import numpy as np
import sys
from pathlib import Path

def readData(jsonLocation, jsonName):
    with open(os.path.join(jsonLocation, jsonName), 'r') as file:
        data = json.load(file)
    return data

def saveData(jsonLocation, jsonName, jsonData):
    with open(os.path.join(jsonLocation, jsonName), 'w') as file:
        json.dump(jsonData, file, indent=1)

def removeOvergrownSims(jsonLocation, jsonName, faceFolderLocation, overgrownThresh):
    jsonFile = readData(jsonLocation, jsonName)
    dirList = os.listdir(faceFolderLocation)
    keysToRemove = []

    for ensemble in jsonFile:
        for folderName in dirList:
            if ensemble in folderName:
                for fileName in os.listdir(os.path.join(faceFolderLocation, folderName)):
                    faceFile = np.loadtxt(os.path.join(faceFolderLocation, folderName, fileName), delimiter=",", dtype=float)
                    numFaces = len(faceFile)
                    if numFaces > overgrownThresh:
                        keysToRemove.append(ensemble)
                        print(f"Found {fileName}. Number of faces = {numFaces}.")
                        break
    for key in keysToRemove:
        jsonFile.pop(key)
        saveData(jsonLocation, jsonName, jsonFile)
        print(f"Removed key {key} from {jsonName}.")

def main():
    # Usage: python removeOvergrownSims.py overgrown_threshold_number "path/to/json_directory" "path/to/face_folder_directory" "json_file_name.json"
    # For the (4,2,2,1) toroid, overgrown_threshold_number=96
    userArgs = sys.argv
    overgrownThresh = int(userArgs[1])
    jsonLocation = Path(userArgs[2])
    faceFolderLocation = Path(userArgs[3])
    jsonName = userArgs[4]

    removeOvergrownSims(jsonLocation, jsonName, faceFolderLocation, overgrownThresh)

if __name__=="__main__":
    main()