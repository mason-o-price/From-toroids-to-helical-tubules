import os
import numpy as np
import matplotlib.pyplot as plt
import json
from collections import defaultdict

def SignedVolume(a,b,c,d):
    return 1/6*np.dot(np.cross(b-a,c-a),d-a)

def edgesInFace(face):
    return [[face[0], face[1]], [face[1], face[2]], [face[2], face[0]]]

def edgeFaceIntersect(edge, face, V):
    p1 = V[face[0]]
    p2 = V[face[1]]
    p3 = V[face[2]]

    q1 = V[edge[0]]
    q2 = V[edge[1]]

    if SignedVolume(q1, p1, p2, p3)*SignedVolume(q2, p1, p2, p3)<0 and SignedVolume(q1, q2, p1, p2)*SignedVolume(q1,q2,p2,p3)>0 and SignedVolume(q1,q2,p2,p3)*SignedVolume(q1,q2,p3,p1)>0:
        return True

    # See: https://stackoverflow.com/questions/42740765/intersection-between-line-and-triangle-in-3d
    # If SignedVolume(q1,p1,p2,p3) and SignedVolume(q2,p1,p2,p3) have different signs
    # AND SignedVolume(q1,q2,p1,p2), SignedVolume(q1,q2,p2,p3) and SignedVolume(q1,q2,p3,p1)
    # have the same sign, then there is an intersection.

def testEdgesCrossFace(face1Edges, face2, V):
    for edge in face1Edges:
        if edgeFaceIntersect(edge, face2, V):
            return True
    return False

def faceFaceIntersect(face1, face2, V):
    face1Edges = edgesInFace(face1)
    face2Edges = edgesInFace(face2)

    if testEdgesCrossFace(face1Edges, face2, V):
        return True
    if testEdgesCrossFace(face2Edges, face1, V):
        return True
    return False

def saveCountDataOverTime(countData, outputLocation, outputName):
    with open(os.path.join(outputLocation, outputName), 'w') as f:
        json.dump(countData, f)

def countIntersections(V_location, F_location):
    countData = defaultdict(list)

    timeCharStart = 7
    timeCharEnd = -16
    
    dirList = os.listdir(V_location)
    for folderName in dirList:
        simName = folderName.split("ensemble")[1]
        fileList = os.listdir(os.path.join(V_location, folderName))
        maxIntersections = 0
        for VfileName in fileList:
            if "coordinates" in VfileName:
                timeStep = VfileName[timeCharStart:timeCharEnd]
                V = np.loadtxt(os.path.join(V_location, folderName, VfileName), delimiter=",", dtype=float)
                faceFilename = VfileName.replace("lammps", "vtk")
                faceFilename = faceFilename.replace("coordinates", "faces")
                F = np.loadtxt(os.path.join(F_location, folderName, faceFilename), delimiter=",").astype(np.int64)

                count = 0

                for (idx, face1) in enumerate(F):
                    if isinstance(face1, np.ndarray) :
                        for face2 in F[idx+1:]:
                            if faceFaceIntersect(face1, face2, V):
                                count += 1
                    else:
                        break
                countData[simName].append((timeStep, count))
                if count > maxIntersections:
                    maxIntersections = count
        print(f"Processed {folderName}. Found a max of {maxIntersections} intersections. Finished with {count} intersections.")
    return countData

def main():
    V_location = r'C:\Users\mason\OneDrive\Desktop\code\openMesh\test_20250928_toroid\data_4221_toroid_on_target_R_exc0_165_med_res_coordinates' # vertex coordinates file
    F_location = r'C:\Users\mason\OneDrive\Desktop\code\openMesh\test_20250928_toroid\data_4221_toroid_on_target_R_exc0_165_med_res_faces' # faces

    outputLocation = r'C:\Users\mason\OneDrive\Desktop\code\openMesh\analysis\intersectionAnalysis\intersectionData'
    outputName = "intersections_4221_toroid_R_exc0_165_on_target_med_res.json"

    countData = countIntersections(V_location, F_location)
    saveCountDataOverTime(countData, outputLocation, outputName)

if __name__ == "__main__":
    main()
    