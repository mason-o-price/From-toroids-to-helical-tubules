import matplotlib.pyplot as plt
import json
import numpy as np
import os
import sys
from pathlib import Path

def readData(dataLocation, fileName):
    with open(os.path.join(dataLocation, fileName), 'r') as file:
        data = json.load(file)
    return data

def unpackCountData(sim, countData):
    Xdata = [float(t) for t,_ in countData[sim]]
    Ydata = [count for _,count in countData[sim]]
    zipped = sorted(zip(Xdata, Ydata), key=lambda x: x[0])
    Xdata, Ydata = zip(*zipped)
    Xdata, Ydata = list(Xdata), list(Ydata)
    return Xdata, Ydata

def plotCountsOverTime(countData, yMin, yMax, lineWidth, lineAlpha):
    for sim in countData:
        Xdata, Ydata = unpackCountData(sim, countData)
        plt.plot(Xdata, Ydata, color='blue', linewidth=lineWidth, alpha=lineAlpha)
    plt.ylim(yMin, yMax)
    plt.show()

def countFramesWithIntersections(countData):
    intersectionRates = []
    for sim in countData:
        countList = np.array([count for _,count in countData[sim]])
        totalFrames = len(countList)
        intersectingFrames = np.count_nonzero(countList)
        intersectionRates.append(intersectingFrames/totalFrames)
    return np.average(np.array(intersectionRates))

def findAverageIntersectionsPerSim(countData):
    averageIntersections = 0
    numSimulations = len(countData)
    for sim in countData:
        _, Ydata = unpackCountData(sim, countData)
        averageIntersections += Ydata[-1]
    averageIntersections /= numSimulations
    return averageIntersections

def main():
    # Usage: python plotCountsOverTime.py "path/to/directory" "fileName.json"
    userArgs = sys.argv
    dataLocation = Path(userArgs[1])
    fileName = userArgs[2]
    simulationData = readData(dataLocation, fileName)
    
    yMax = 20
    yMin = -1
    lineWidth = 1
    lineAlpha = 0.2
    
    plotCountsOverTime(simulationData, yMin, yMax, lineWidth, lineAlpha)
    numFramesWithIntersections = countFramesWithIntersections(countData=simulationData)
    averageIntersectionsPerSim = findAverageIntersectionsPerSim(countData=simulationData)
    print(f"Rates of intersecting frames: {numFramesWithIntersections}")
    print(f"Average intersections per simulation at the final snapshot: {averageIntersectionsPerSim}")

if __name__=="__main__":
    main()