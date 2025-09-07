function uniqueTriangles = findUniqueTriangles(vertexConnectivity, coordinates, numTriangles, m)
%FINDUNIQUETRIANGLES Find the set of unique triangles from vertex connections
%
% uniqueTriangles = findUniqueTriangles(vertexConnectivity, coordinates, numTriangles, m)
%
% INPUT: 
%  vertexConnectivity # cell matrix of connected vertices
%  coordinates # matrix of coordinates
%  numTriangles # number of triangles
%  m # tubelet circumference (T)
% 
% OUTPUT:
%  uniqueTriangles # Cell matrix of unique triangles

trianglesMatrix = cell(numTriangles, 1);
counter = 1;
for k = 1:size(coordinates, 1) %go through all of the points
    point = coordinates(k,1:2);
    point = [mod(point(1)-1, m)+1, point(2)];
    i = point(1); %Find the i,h coordinates
    h = point(2);

    if ~isempty(vertexConnectivity{i,h}) %check if the point is deleted -> if the point's cell is not empty, continue...
        %Find the neighbors
        neighbors = vertexConnectivity{i,h};
        for j = 1:size(neighbors,1)
            n = neighbors(j, :); %This will give you the i,h values for this neighbor

            %If two vertices share another vertex in common --> they form a triangle
            n_neighbors = vertexConnectivity{n(1), n(2)}; %find the other point's neighbors
            intersection = intersect(n_neighbors, neighbors, 'rows'); %Find the intersection of the cells. This could be more than one vector.
            %intersection = intersection(1,:); % take only the first row in the intersection
            if ~isempty(intersection) %If the two vertices share a neighbor in common, draw a patch of a plane bounded by the 3 vertices
                for row = 1:size(intersection,1)
                    final_vertex = intersection(row, :);
                    if ~isequal([h, h, h], [h, n(2), final_vertex(2)]) %Make sure that we aren't filling the cross-sections of a tubule, (only happens for a (3,0) tubule).
                        trianglesMatrix{counter, 1} = [i,h; n(1), n(2); final_vertex(1), final_vertex(2)];
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end

% Remove redundant triangles
uniqueTriangles = cell(numTriangles, 1);
uniqueTriangles{1,1} = trianglesMatrix{1,1}; %copy the first cell
counter = 1;
for row = 1:size(trianglesMatrix, 1)
    add_cell = 1; % truth value of whether we should add the cell to the unique triangles matrix
    cell_k = trianglesMatrix{row, 1}; %this is the cell to compare
    for j = 1:size(uniqueTriangles, 1)
        cell_j = uniqueTriangles{j, 1}; %compare it to all of the cells in unique triangles matrix
        if isequal(sort(cell_k), sort(cell_j))
            add_cell = 0;
        end
    end
    if add_cell
        counter = counter + 1;
        uniqueTriangles{counter, 1} = cell_k;
    end
end

%% Find the orientations of the triangles & reorder the matrices
for row = 1:size(uniqueTriangles,1)
    triangle = uniqueTriangles{row, 1};
    if ~isempty(triangle)
        for j = 1:3
            vertex = triangle(j, :);
            if vertex(2) ~= mode(triangle(:,2)) % if the vertex is the top-most of an upward triangle, or the bottom-most of a downward triangle
                uniqueTriangles{row,1} = circshift(triangle, 4-j); %then move that point to the top of the triangle's matrix
            end
        end
    end
end


end