function new_triangles_matrix = reorderTriangles(triangleMatrix, T,L,D,R,N)
%REORDERTRIANGLES Simplify the order of the triangles matrix
%
% new_triangles_matrix = REORDERTRIANGLES(triangleMatrix, c,L,d,x,N)
%
% INPUT:
%  triangleMatrix # Initial triangles matrix
%  T # Tubelet circumference
%  L # Tubelet length
%  D # Junction side-length
%  R # Rolling parameter
%  N # Number of tubelet segments
%
% OUTPUT:
%  new_triangles_matrix # Update matrix of triangles

% Update the name of the variables
c = T;
d = D;
x = R;

new_triangles_matrix = cell(size(triangleMatrix, 1), 1);
counter = 1;
for j = 1:N*(L+d)
    for i = 1:2*c
        if mod(i,2) ~= 0 % triangle pointing upward
            marker_vtx = [mod(i-1+x*floor((j-1)/(L+d))-floor((i-1)/2), c)+1, j+1];
            for k = 1:size(triangleMatrix, 1)
                triangle = triangleMatrix{k, 1};
                if isequal(triangle(1,:), marker_vtx) && triangle(2,2) < triangle(1,2) % check that triangle includes the marker vertex, and that it is pointing up
                    new_triangles_matrix{counter, 1} = triangle;
                    counter = counter + 1;
                end
            end
        elseif mod(i,2) == 0 % triangle pointing downward
            marker_vtx = [mod(i-1+x*floor((j-1)/(L+d))-floor((i-1)/2), c)+1, j];
            for k = 1:size(triangleMatrix, 1)
                triangle = triangleMatrix{k,1};
                if isequal(triangle(1,:), marker_vtx) && triangle(2,2) > triangle(1,2)
                    new_triangles_matrix{counter, 1} = triangle;
                    counter = counter + 1;
                end
            end
        end
    end
end

end