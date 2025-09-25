import matplotlib.pyplot as plt
import json
import numpy as np
import os

def readData(dataLocation, fileName):
    with open(os.path.join(dataLocation, fileName), 'r') as file:
        data = json.load(file)
    return data

def plotCountsOverTime(countData, yMin, yMax, lineWidth, lineAlpha):
    plt.figure()
    for sim in countData:
        Xdata = [float(t) for t,_ in countData[sim]]
        Ydata = [count for _,count in countData[sim]]
        zipped = sorted(zip(Xdata, Ydata), key=lambda x: x[0])
        Xdata, Ydata = zip(*zipped)

        Xdata, Ydata = list(Xdata), list(Ydata)
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
    return intersectionRates

def main():
    dataLocation = r"C:\Users\mason\OneDrive\Desktop\code\openMesh\analysis\intersectionAnalysis\intersectionData"
    fileName = "intersections_4221_toroid_R_exc0_15_on_target_high_res.json"
    simulationData = readData(dataLocation, fileName)

    yMax = 10
    yMin = -1
    lineWidth = 1
    lineAlpha = 0.25
    
    plotCountsOverTime(simulationData, yMin, yMax, lineWidth, lineAlpha)
    numFramesWithIntersections = countFramesWithIntersections(countData=simulationData)
    print(f"Rates of intersecting frames: {numFramesWithIntersections}")



if __name__=="__main__":
    main()