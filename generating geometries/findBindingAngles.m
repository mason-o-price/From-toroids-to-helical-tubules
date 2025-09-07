function [cell_matrix, angles] = findBindingAngles(coordinates, trianglesMatrix, speciesVector, numSpecies, m)
%FINDBINDINGANGLES Measure the dihedral angle between neighboring triangles
%
% [cell_matrix, angles] = FINDBINDINGANGLES(coordinates, trianglesMatrix, speciesVector, numSpecies, m)
%
% INPUT:
%  coordinates # vertex coordinates
%  trianglesMatrix # Cell matrix storing the triangles. Entry at {i,j} = [v1; v2; v3] 
%  speciesVector # Nx1 vector of species for each triangle
%  numSpecies # integer number of species
%  m # tubule circumference
%
% OUTPUT:
% cell_matrix # cell matrix storing the binding angles between all instances of different species
% angles # matrix of mean binding angles between species

cell_matrix = cell(numSpecies, numSpecies);
angles = zeros(numSpecies, numSpecies);
% format: row = species of triangle 1; column = species of triangle 2.

for j = 1:size(trianglesMatrix, 1)
    triangle_A = trianglesMatrix{j,1}; % find the first triangle
    for k = 1:size(trianglesMatrix, 1)
        if k ~=j
            triangle_B = trianglesMatrix{k,1}; % find the second triangle
            if size(intersect(triangle_A,triangle_B, 'rows'), 1) == 2 % if they share 2 vertices in common --> they are neighbors
                % Determine the vertices (& their order) in each triangle.
                % Do triangle A:
                v1_A = triangle_A(1,:); % find vertex 1
                v2_A = triangle_A(2,:); % vertex 2
                v3_A = triangle_A(3,:); % vertex 3

                % Now do triangle B:
                v1_B = triangle_B(1,:); % find vertex 1
                v2_B = triangle_B(2,:); % vertex 2
                v3_B = triangle_B(3,:); % vertex 3
             
                v1_A_XYZ = coordinates(findIndx(v1_A(1), v1_A(2), m), 3:5); %find the xyz coordinates of v1_A
                v2_A_XYZ = coordinates(findIndx(v2_A(1), v2_A(2), m), 3:5); %find the xyz coordinates of v2_A
                v3_A_XYZ = coordinates(findIndx(v3_A(1), v3_A(2), m), 3:5); %find the xyz coordinates of v3_A
                v1_B_XYZ = coordinates(findIndx(v1_B(1), v1_B(2), m), 3:5); %find the xyz coordinates of v1_B
                v2_B_XYZ = coordinates(findIndx(v2_B(1), v2_B(2), m), 3:5); %find the xyz coordinates of v2_B
                v3_B_XYZ = coordinates(findIndx(v3_B(1), v3_B(2), m), 3:5); %find the xyz coordinates of v3_B
                
                % Determine the (counter-clockewise) oriented edge vectors on triangle A
                e1_A = v1_A_XYZ - v3_A_XYZ;
                e2_A = v2_A_XYZ - v1_A_XYZ;
                e3_A = v3_A_XYZ - v2_A_XYZ;

                % Determine the (counter-clockwise) oriented edge vectors on triangle B
                e1_B = v1_B_XYZ - v3_B_XYZ;
                e2_B = v2_B_XYZ - v1_B_XYZ;
                e3_B = v3_B_XYZ - v2_B_XYZ;

                % Determine the unit face normal vectors:
                N_A = cross(e1_A, e2_A)/norm(cross(e1_A, e2_A));
                N_B = cross(e1_B, e2_B)/norm(cross(e1_B, e2_B));

                center_A = mean([v1_A_XYZ; v2_A_XYZ; v3_A_XYZ]); % Find the center of triangle A
                center_B = mean([v1_B_XYZ; v2_B_XYZ; v3_B_XYZ]); % Find the center of triangle B
                
                d1 = norm(center_A - center_B);
                
                offset_A = center_A + 0.2*N_A; % follow the normal vector away from the center of triangle A
                offset_B = center_B + 0.2*N_B; % follow the normal vector away from the center of triangle B

%                % If checking offset locations
%                plot3(offset_A(1), offset_A(2), offset_A(3), 'b*');
%                plot3(offset_B(1), offset_B(2), offset_B(3), 'ro');

                d2 = norm(offset_A - offset_B);

                if d2 < d1 % if the angle is negative (concave), then the offset vectors should become closer.
                    sign = -1; % set the angle to be negative.
                else %otherwise, assume the angle is positive.
                    sign = 1;
                end

                % Calcualte the binding angle: theta = atan2(
                bindingAngle = 180/pi*sign*abs(atan2(norm(cross(N_A, N_B)), dot(N_A, N_B)));

                row = speciesVector(j); % species of triangle 1 = row in angles matrix
                column = speciesVector(k); % species of triangle 2 = column in angles matrix
    
                cell_matrix{row, column} = [cell_matrix{row, column}; bindingAngle];

            end
        end
    end

for i = 1:size(cell_matrix,1)
    for j = 1:size(cell_matrix,2)
        possible_angles = cell_matrix{i,j};
        if isempty(possible_angles)
            angles(i,j) = NaN;
        else
            angles(i,j) = mean(possible_angles);
        end
    end
end

end
