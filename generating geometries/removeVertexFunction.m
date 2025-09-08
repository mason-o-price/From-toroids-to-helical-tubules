function vertexConnectivity = removeVertexFunction(i,h,vertexConnectivity)
%REMOVEVERTEXFUNCTION Delete vertex and remove from neighbors' connections
%
% vertexConnectivity = REMOVEVERTEXFUNCTION(i,h,vertexConnectivity)
%
% INPUT:
%  i # lattice coordinate 1
%  h # lattice coordinate 2
%  vertexConnectivity # Initial connections between vertices. 
%
% OUTPUT:
%  vertexConnectivity # Updated connections between vertices.


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