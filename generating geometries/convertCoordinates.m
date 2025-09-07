function [xy_matrix] = convertCoordinates(hk_matrix)
%CONVERTCOORDINATES Convert lattice vectors (h,k) into cartesion coordinates (x,y)
%
% xy_matrix = CONVERTCOORDINATES(hk_matrix)
%
% INPUT:
%  hk_matrix # Nx2 matrix of (h,k) values for N points
%
% OUTPUT:
%  xy_matrix # Nx2 matrix of (x,y) coordinates

xy_matrix = [];
for row = 1:size(hk_matrix, 1)
    h = hk_matrix(row, 1);
    k = hk_matrix(row, 2);
    x = h + k/2;
    y = k*sqrt(3)/2;
    xy_matrix(row, 1) = x;
    xy_matrix(row, 2) = y;
end