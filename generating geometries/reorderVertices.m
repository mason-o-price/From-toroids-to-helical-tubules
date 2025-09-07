function newTrianglesMatrix = reorderVertices(triangleMatrix, T)
%REORDERVERTICES Simplify the order of the vertices in the triangles matrix
%
% newTrianglesMatrix = REORDERVERTICES(triangleMatrix, c)
%
% INPUT:
%  triangleMatrix # Matrix of triangles
%  T % Tubelet circumference
%
% OUTPUT:
%  newTrianglesMatrix # Updated list of triangles with re-ordered vertices

newTrianglesMatrix = triangleMatrix;

for k = 1:size(triangleMatrix, 1)
    triangle = triangleMatrix{k, 1}; % go through each triangle
    i2 = triangle(2,1); % find the horizontal lattice coordinates of the second and third vertices in the matrix
    i3 = triangle(3,1); 
    if triangle(2,2) < triangle(1,2)% if the triangle is pointing up
        % handle the case along the seem of the tubule, vertex c should come before 1
        if ( i2 == 1 ) && ( i3 == T ) % in this case, swap the order of v2 & v3, which would not otherwise happen in the normal case
            triangle([2,3], :) = triangle([3,2], :);
            newTrianglesMatrix{k, 1} = triangle;
        elseif i2 > i3  % otherwise, treat it like a normal vertex & swap v2 & v3 only if v3 is ordered before v2
            triangle([2,3], :) = triangle([3,2], :);
            newTrianglesMatrix{k, 1} = triangle;
        end
    end
    if triangle(2,2) > triangle(1,2) % if the triangles is pointing downif v3 is ordered before v2 (i.e. vertex on the left comes before the vertex on the right), swap the rows
        % Handle the case along the seem of the tubule, vertex 1 should come before c.
        % Ordinarily, the vertices would be swapped under this condition. Avoid swapping if this condition is met.
        if ( i2 < i3 ) && ~( ( i2 == 1 ) && ( i3 == T ) ) % If they appear out of order, and we're not along a seem, swap v2 & v3
            triangle([2,3], :) = triangle([3,2], :);
            newTrianglesMatrix{k, 1} = triangle;
        end
    end
end

end