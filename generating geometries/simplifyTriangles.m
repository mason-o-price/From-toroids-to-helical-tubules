function trianglesNoCell = simplifyTriangles(trianglesMatrix, T)
%SIMPLIFYTRIANGLES Convert a cell matrix of triangles into a regular matrix
%
% trianglesNoCell = SIMPLIFYTRIANGLES(trianglesMatrix, T)
%
% INPUT:
%  trianglesMatrix # Cell-matrix of triangles
%  T # Tubelet circumference
%
% OUPUT:
%  trianglesNoCell # Strandard matrix of triangle faces, row = [v1,v2,v3]

numTriangles = size(trianglesMatrix, 1);
trianglesNoCell = zeros(numTriangles, 3);

for j = 1:numTriangles
    triangle = trianglesMatrix{j,1};
    v1 = triangle(1,:);
    v2 = triangle(2,:);
    v3 = triangle(3,:);

    v1Indx = findIndx(v1(1), v1(2), T);
    v2Indx = findIndx(v2(1), v2(2), T);
    v3Indx = findIndx(v3(1), v3(2), T);

    trianglesNoCell(j,:) = [v1Indx, v2Indx, v3Indx];
end
