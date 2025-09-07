function plotFaces(coordinates, trianglesMatrix, m, plot_markers, speciesVector, colorMatrix) 
%PLOTFACES Plot the set of edges for a curved tubule.
% 
% PLOTFACES(coordinates, vertexConnections, m, edgeColor)
%
% INPUT:
%  coordinates # Matrix of vertex coordinates
%  trianglesMatrix # Matrix of triangles
%  m # Tubule circumference (T)
%  plot_markers # whether to plot edge-orientation markers
%  speciesVector # vector of triangle species
%  colorMatrix # Matrix of face colors
%
% NOTES:
%  If no colors/species --> use 0 for each parameter.

%In case if we don't want to plot colors:
if speciesVector == 0
    speciesVector = ones(size(coordinates, 1));
end
if colorMatrix == 0
    colorMatrix = ones(size(coordinates, 1), 3);
end

for row = 1:size(trianglesMatrix, 1) 
    triangle = trianglesMatrix{row, 1}; % each row corresponds to a triangle
    species = speciesVector(row);
    if species == 0
        color = [1,1,1];
    else
        color = colorMatrix(species, :);
    end
    XYZ = [];
    for j = 1:3
        XYZ = [XYZ; coordinates(findIndx(triangle(j,1), triangle(j,2), m), 3:5)];
    end
    %species = speciesMatrix(row,1);
    patch(XYZ(:,1), XYZ(:,2), XYZ(:,3), color);%, colorMatrix(species,:));

    %plot the marker line-segments
    if plot_markers
        center_pnt = mean(XYZ);
        lines = line([center_pnt(1); XYZ(1,1)],...
            [center_pnt(2); XYZ(1,2)],...
            [center_pnt(3); XYZ(1,3)]); 
        set(lines, 'color', 'k');
    end
end
end