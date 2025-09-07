function connectionsNoCell = simplifyConnections(vertexConnectivity, coordinates, T)
%SIMPLIFYCONNECTIONS Convert a cell matrix of vertex connections into a
%regular matrix
% 
% connectionsNoCell = SIMPLIFYCONNECTIONS(vertexConnectivity, coordinates, T)
%
% INPUT:
%  vertexConnectivity # Cell-matrix of vertex connections
%  coordinates # Matrix of vertex coordinates
%  T # Tubelet circumference
%
% OUPUT:
%  connectionsNoCell # Regular matrix of vertex connections (index = row in
%  coordinates matrix)

numVertices = size(coordinates, 1);
connectionsNoCell = zeros(numVertices, 8);

for j = 1:numVertices
    vertex = coordinates(j,1:2);
    neighbors = vertexConnectivity{vertex(1), vertex(2)};
    for k = 1:size(neighbors, 1)
        vk = neighbors(k,:);
        vkIndx = findIndx(vk(1), vk(2), T);
        connectionsNoCell(j,k) = vkIndx;
    end
end
