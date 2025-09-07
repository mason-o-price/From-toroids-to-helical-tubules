function [coordinates, vertexConnectivity, trianglesMatrix, kMatrix] = mergeHoles(coordinates, vertexConnectivity, trianglesMatrix, Pcollection, Qcollection, m, kMatrix, kConstTubule)
%MERGEHOLES Merge the holes across a holey tubule
% 
% [coordinates, vertexConnectivity, trianglesMatrix, kMatrix] = mergeHoles(coordinates, vertexConnectivity, trianglesMatrix, Pcollection, Qcollection, m, kMatrix, kConstTubule)
% 
% INPUT:
%  coordinates # Matrix of vertex coordinates
%  vertexConnectivity # Cell matrix of connected vertices
%  trianglesMatrix # Matrix of triangles
%  Pcollection # Set of starting vertices in pairs to be merged together
%  Qcollection # Set of terminal vertices in pairs to be merged together
%  m # Tubelet circumference (T)
%  kMatrix # Matrix of spring constants
%  kConstTubule # Constant spring force for edges within the tubule
%
% OUTPUT:
%  coordinates # Updated set of coordinates
%  vertexConnectivity # Updated set of vertex connections
%  trianglesMatrix # Updated matrix of triangles
%  kMatrix # Updated matrix of spring constants

% Stitch the vertices across the holes (to update the neighbors)
for row = 1:size(Pcollection, 1)
    p = Pcollection(row, :); %find the first point
    q = Qcollection(row, :); %find the second point
    
    % Merge the vertices
    [coordinates, vertexConnectivity, kMatrix] = mergeVertex(p, q, coordinates, vertexConnectivity, Qcollection, m, kMatrix, kConstTubule);
end

% Delete all of the vertices that we already stitched up & replace them
for i = 1:size(Qcollection, 1)
    q = Qcollection(i,:);
    vertexConnectivity = removeVertex(q(1), q(2), vertexConnectivity);
    for row = 1:size(trianglesMatrix,1) % redefine the triangles in the triangles matrix
        triangle = trianglesMatrix{row, 1}; 
        [is_a_vertex, vertex_location] = ismember(q, triangle, 'rows'); % check if this point is in the triangle, and if so find its location
        if is_a_vertex
            p = Pcollection(i,:); 
            triangle(vertex_location,:) = p; %replace the deleted point with it's associated vertex
            trianglesMatrix{row,1} = triangle; % update the triangles matrix
        end
    end
    q_indx1 = (q(2)-1)*m+q(1);
    kMatrix(:,q_indx1) = zeros(size(kMatrix, 1), 1);
    kMatrix(q_indx1,:) = zeros(1, size(kMatrix, 2)); %give any left-over bonds a strength of zero

end


%% Define the mergeVertex function
function [coordinates, vertexConnectivity, kMatrix] = mergeVertex(p,q,coordinates, vertexConnectivity, Qcollection, m, kMatrix, kConstTubule)
%Update the coordinates and connectivity matrix after merging points p,q

p_cell = vertexConnectivity{p(1), p(2)}; %Grab the cell of the first vertex
q_cell = vertexConnectivity{q(1), q(2)}; %Grab the cell of the second vertex

% Form new bonds
if ~isempty(p_cell) && ~isempty(q_cell)
    for row = 1:size(q_cell,1) %Go through the cell of q
        neighbor = q_cell(row, :); %Consider each neighbor
        %Check if the neighbor is in the set of vertices that we need to delete (Qcollection)
        shouldConnect = 1; %Should we connect the two vertices? Initially say yes.
        for j = 1:size(Qcollection,1)  
            if neighbor == Qcollection(j,:)
                shouldConnect = 0; %If the proposed point is in the set of vertices to delete, don't connect.
            end
        end
        for j = 1:size(p_cell,1) %check if we're connecting p to a point that is already its neighbor.
            if neighbor == p_cell(j,:)
                shouldConnect = 0;
            end
        end
        if shouldConnect
            %update the vertex connecitivity matrix
            vertexConnectivity{p(1), p(2)} = [vertexConnectivity{p(1), p(2)}; neighbor]; %Add the neighbor to the list of neighbors of p
            vertexConnectivity{neighbor(1), neighbor(2)} = [vertexConnectivity{neighbor(1), neighbor(2)}; p]; % connect both ways.
        
            %update the kMatrix storing the spring constants for edges
            kMatrix(findIndx(p(1), p(2), m), findIndx(neighbor(1), neighbor(2), m)) = kConstTubule;
            kMatrix(findIndx(neighbor(1), neighbor(2), m), findIndx(p(1), p(2), m)) = kConstTubule;
        end
    end

    p_indx = findIndx(p(1), p(2), m); %Find the index of the points in the coordinates matrix
    q_indx = findIndx(q(1), q(2), m);

    p_coordinates = coordinates(p_indx, 3:5); %Grab the xyz coordinates of the points
    q_coordinates = coordinates(q_indx, 3:5); 

    avgPosition = (p_coordinates + q_coordinates)/2; %Take the average of their positions

    coordinates(p_indx, 3:5) = avgPosition; %Redefine the position of p to be the average of p & q
end
end

%% Define the remove vertex function
function [vertexConnectivity] = removeVertex(i,h,vertexConnectivity)
%This will delete the vertex at (i,h) from the vertexConnectivity &
%uniqueConnections matrices, it will also remove them from their neighbor's list of neighors.

%Go into the cell associated with this vertex
vertex_cell = vertexConnectivity{i, h};
if ~isempty(vertex_cell)
    for j = 1:size(vertex_cell,1) %go through each of its neighbors
        neighbor = vertex_cell(j,1:2); %find the values of i & h of the jth neighbor
        
        %We need to remove our vertex from each of its neighbor's list of neighbors
        neighbor_cell = vertexConnectivity{neighbor(1), neighbor(2)}; %Access the cell of this neighbor
        if ~isempty(neighbor_cell)
            length = size(neighbor_cell(:,1));
            %check until we find our vertex, and remove it.
            for row = 0:length(1)-1 %we need to check the vector in reverse order because we might remove stuff before we get to the end, which will mess up the indexing.
                if neighbor_cell(length(1) - row,:) == [i, h]
                    vertexConnectivity{neighbor(1), neighbor(2)}(length(1) - row,:) = []; %remove it from vertexConnectivity
                end
            end
        end
    end
end

%Remove the cell of the vertex that we want to delete
vertexConnectivity{i, h} = [];
end

end