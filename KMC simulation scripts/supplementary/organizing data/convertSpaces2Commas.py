import pandas as pd
import os
import sys
from pathlib import Path

def convertSpaces2Commas(inputDirectory, outputDirectory):
    for fileName in os.listdir(inputDirectory):
        print(f"Processing {fileName}")
        df = pd.read_csv(os.path.join(inputDirectory, fileName), sep=' ', header=None) # Read file
        df = df.iloc[:, 1:] # remove first column
        df.to_csv(os.path.join(outputDirectory, fileName), sep=',', index=False, header=False)

def main():
    # Usage: python convertSpaces2Commas.py "path/to/input_directory" "path/to/output_directory"
    userArgs = sys.argv
    inputDirectory = Path(userArgs[1])
    outputDirectory = Path(userArgs[2])

    convertSpaces2Commas(inputDirectory,outputDirectory)

if __name__=="__main__":
    main()