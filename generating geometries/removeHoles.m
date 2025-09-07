function vertexConnectivity = removeHoles(m, x, c, numTubules, tubuleHeight, junctionSideLength, vertexConnectivity)
%REMOVEHOLES Remove holes of specified size and location in a tubule
%
% vertexConnectivity = REMOVEHOLES(m, x, c, numTubules, tubuleHeight, junctionSideLength, vertexConnectivity)
%
% INPUT:
%  m # Tubule circumference (T)
%  x # Spacing between cuts (R)
%  c # Constant offset in h-lattice direction
%  numTubules # Number of tubelet segments (N)
%  tubuleHeight # Total number of triangles along the height of the tubule
%  junctionSideLength # Side-length of junction region (D)
%  vertexConnectivity # Matrix of vertex connections
%
% OUTPUT:
%  vertexConnectivity # Updated matrix of vertex connections

%Define a function to remove an individual vertex
removeVertex = @(i,h,vertexConnectivity)removeVertexFunction(i,h,vertexConnectivity);
%makePeriodic = @(vertex)makePeriodic(vertex, m);

%Remove all of the vertices in the interior of a hole.
for k = 1:numTubules
    %Define the coordinates of the bottom right corner of the junction
    a2 = [c + (k-1)*x + junctionSideLength, k*tubuleHeight + (k-1)*junctionSideLength];
    a2 = makePeriodic(a2, m);
    a2(2) = a2(2) + 1; %Offset the height by 1

    %Go through all the vertices that need to be deleted
    for u = 1:m-1-junctionSideLength
        for v = 1:junctionSideLength-1
            point_to_remove = makePeriodic(a2 + [u,v], m); %the i value of the vertex to be deleted
            %remove the vertex; update the connections matrices. 
            vertexConnectivity = removeVertex(point_to_remove(1), point_to_remove(2), vertexConnectivity);
        end
    end
end

%% Remove the corners of the hole (around the NGV defect points)
p_matrix = [];
q_matrix = []; %matrices to store the vertex pairs
for k = 1:numTubules
    %Define the coordinates of the bottom right corner of the junction
    a2 = makePeriodic([c + (k-1)*x + junctionSideLength, k*tubuleHeight + (k-1)*junctionSideLength], m);
    a2 = a2 + [0,1]; %offset the start to match MATLAB's indexing

    %Do the same for the top-left corner
    a4 = makePeriodic([c + (k-1)*x, k*tubuleHeight + k*junctionSideLength], m);
    a4 = a4 + [0,1]; %offset the start to match MATLAB's indexing

    %Now find the points across which the corner edges remain
    p1 = makePeriodic(a2 + [1,0], m); 
    q1 = makePeriodic(a2 + [0,1], m);
    p2 = makePeriodic(a4 - [1, 0], m); 
    q2 = makePeriodic(a4 - [0, 1], m);
    %p1 <-> q1
    %p2 <-> q2
    p_matrix = [p_matrix; p1];
    p_matrix = [p_matrix; p2];
    q_matrix = [q_matrix; q1];
    q_matrix = [q_matrix; q2];

    %If junctionSideLength is equal to 1, we need to remove additional
    %extraneous edges.
    if junctionSideLength == 1
        for i = 1:m-2
            p = a2 + [i,0];
            p = makePeriodic(p, m);
            q = p + [0, junctionSideLength];
            q = makePeriodic(q, m);
            p_matrix = [p_matrix; p];
            q_matrix = [q_matrix; q];
        end
    end
end
%Now remove them from each other's list of neighbors.
%Remove p1 from q1's neighbors
for indx = 1:size(p_matrix,1)
    p = p_matrix(indx, :);
    q = q_matrix(indx, :);
    for row = 1:size(vertexConnectivity{q(1), q(2)},1)
        if vertexConnectivity{q(1), q(2)}(row, :) == p
            vertexConnectivity{q(1), q(2)}(row, :) = [];
            break
        end
    end

    %Remove q1 from p1's neighbors
    for row = 1:size(vertexConnectivity{p(1), p(2)},1)
        if vertexConnectivity{p(1), p(2)}(row, :) == q
            vertexConnectivity{p(1), p(2)}(row, :) = [];
            break
        end
    end
end

%% Remove the extra ring of vertices from the very bottom
%Define the position of the starting corner of the hole
a3_last = [c + (numTubules-1)*x + junctionSideLength, numTubules*tubuleHeight + numTubules*junctionSideLength];
a3_last = makePeriodic(a3_last, m); %implement the periodic boundary
a3_last = a3_last + [0,1]; %offset the height, I forget why we need to do this.

for i = 1:m-junctionSideLength-1
    point_to_remove = a3_last + [i,0];
    point_to_remove = makePeriodic(point_to_remove, m);
    vertexConnectivity = removeVertex(point_to_remove(1), point_to_remove(2), vertexConnectivity);
end

%% Plot the tubule with the holes removed now


%plotEdges(coordinates, vertexConnectivity, m, 'k');

%% Define the function to remove the vertices

function vertexConnectivity = removeVertexFunction(i,h,vertexConnectivity)
%This will delete the vertex at (i,h) from the vertexConnectivity &
%uniqueConnections matrices, it will also remove them from their neighbor's
%list of neighors.

%Go into the cell associated with this vertex
vertex_cell = vertexConnectivity{i, h};
for j = 1:size(vertex_cell,1) %go through each of its neighbors
    neighbor = vertex_cell(j,1:2); %find the values of i & h of the jth neighbor
    
    %We need to remove our vertex from each of its neighbor's list of neighbors
    neighbor_cell = vertexConnectivity{neighbor(1), neighbor(2)}; %Access the cell of this neighbor
    length = size(neighbor_cell(:,1));
    %check until we find our vertex, and remove it.
    for row = 0:length(1)-1 %we need to check the vector in reverse order because we might remove stuff before we get to the end, which will mess up the indexing.
        if neighbor_cell(length(1) - row,:) == [i, h]
            vertexConnectivity{neighbor(1), neighbor(2)}(length(1) - row,:) = []; %remove it from vertexConnectivity
        end
    end
end

%Remove the cell of the vertex that we want to delete
vertexConnectivity{i, h} = [];

end
end