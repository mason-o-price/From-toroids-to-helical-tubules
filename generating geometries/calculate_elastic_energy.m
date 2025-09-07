function U = calculate_elastic_energy(vertexConnectivity, coordinates, m)
%CALCULATE_ELASTIC_ENERGY Calculate the elastic energy  of a triangular
%mesh for toroidSimulation.m
%
% U = CALCULATE_ELASTIC_ENERGY(vertexConnectivity, coordinates, m)
%
% INPUT:
%  vertexConnectivity # cell matrix storing connected vertices
%  coordinates # vertex coordinates
%  m # tubule circumference
%
% OUTPUT:
%  U # elastic energy assuming the spring constant is normalized.

displacement_matrix = [];

for h = 1:size(vertexConnectivity,1) %go through the rows of the vertex connections matrix
    for k = 1:size(vertexConnectivity,2) %%go through the columns of the vertex connections matrix
        neighbors = vertexConnectivity{h,k}; %find the neighbors
        v_indx = (k-1)*m + h; %find the index of the vertex in the coordinates matrix
        if ~isempty(vertexConnectivity{h,k}) %make sure the first vertex is not deleted
            vertex = coordinates(v_indx, :); %find the coordinates of the vertex. format: cols 1,2,3,4,5 = h,k,x,y,z
            for j = 1:size(neighbors,1)
                q = neighbors(j,:);
                if ~isempty(vertexConnectivity{q(1),q(2)}) %make sure the neighbor is not deleted
                    q_indx = (q(2)-1)*m + q(1); %find the index of the neighbor in the coordinates matrix
                    qXYZ = coordinates(q_indx, 3:5);
                    relativeVector = vertex(3:5) - qXYZ; %find the vector connecting them
                    displacement = relativeVector - relativeVector/norm(relativeVector); %%Find the displacement of the vertices from the rest length (which is the diameter of the tubule here)
                    displacement_matrix = [displacement_matrix; displacement];
                end
            end
            vertexConnectivity = removeVertexFunction(h,k, vertexConnectivity); % delete the first vertex after measuring so that we don't double-count any stretched bonds.
        end
    end
end

displacement_matrix = rmmissing(displacement_matrix); % get rid of all NaN's
U = sum(sqrt(sum(displacement_matrix.^2, 2)),1);