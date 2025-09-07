function speciesVector = findSpecies(T,L,D,R, num_triangles)
%FINDSPECIES Find the species of each triangle in a list of triangles
% 
% speciesVector = findSpecies(c,L,d,x, num_triangles)
% 
% INPUT:
%  T # Tubelet circumference
%  L # Tubelet length
%  D # Junction side-length
%  R # Rolling parameter
%  num_triangles # Number of triangle types
%
% OUTPUT:
%  speciesVector # Vector of species. Entry at j = species of triangle j.


%% Preparation
% Update the variable names
c = T;
d = D;
x = R;

k = c*L+d^2; %number of species

speciesMatrix = zeros(2*c, 2*(L+d)); %tells you which triangles get assigned which colors

offset_x = -1/2;
offset_y = -sqrt(3)/2;

% generate a triangles matrix
trianglesMatrix_2d = cell(c,L+d);
center_of_mass_matrix = cell(size(trianglesMatrix_2d));
for i = 1:2*c
    for h = 1:L+d
        if mod(i,2) ~= 0
            trianglesMatrix_2d{i, h} = [i+(h-1)/2-(i-1)/2, (h+1)*sqrt(3)/2;...
                i+h/2-1-(i-1)/2, h*sqrt(3)/2;...
                i+h/2-(i-1)/2, h*sqrt(3)/2] + [offset_x, offset_y];
            
            center = mean(trianglesMatrix_2d{i,h}, 1);
            center_of_mass_matrix{i,h} = unitPeriodic_xy(center, c,L,d);
        else
            trianglesMatrix_2d{i, h} = [i+h/2-1-i/2+1, h*sqrt(3)/2;...
                i+(h-3)/2-i/2+1, (h+1)*sqrt(3)/2;...
                i+(h-1)/2-i/2+1, (h+1)*sqrt(3)/2] + [offset_x, offset_y];
            center = mean(trianglesMatrix_2d{i,h}, 1);
            center_of_mass_matrix{i,h} = unitPeriodic_xy(center, c,L,d);
        end
    end
end

%Find the symmetry points
symmetry_pnts = [d/2,  L+d/2; ...
    (c+d)/2, L+d/2; ...
    (d-x)/2,  L/2; ...
    (c+d-x)/2, L/2];
symmetry_pnts_xy = convertCoordinates(symmetry_pnts); %+ [offset_x, offset_y];

% assign the right species
counter = 1;
for n = 1:size(center_of_mass_matrix, 2)
    for i = 1:size(center_of_mass_matrix, 1)
        pnt_1 = center_of_mass_matrix{i,n};
        if ( speciesMatrix(i,n) == 0 ) && ( counter <= k ) && (( n <= L ) || ( i <= 2*d )) 
            speciesMatrix(i,n) = counter;
            for j = 1:4
                symmetry_pnt = symmetry_pnts_xy(j,:);
                pnt_2 = unitPeriodic_xy(2*symmetry_pnt-pnt_1, c,L,d,x);
                
                for h_2 = 1:size(trianglesMatrix_2d, 2)
                    for i_2 = 1:size(trianglesMatrix_2d, 1)
                        potential_image = center_of_mass_matrix{i_2, h_2};
                        if norm(potential_image - pnt_2) < 0.01
                            if counter <= k && speciesMatrix(i_2, h_2) ~= counter
                                speciesMatrix(i_2, h_2) = counter;                                    
                            end
                        end
                    end
                end
            end
            counter = counter + 1;
        end
    end
end

%% Cut out holes
punctured_triangles_matrix = [];
for i = 1:size(trianglesMatrix_2d,1)
    for h = 1:size(trianglesMatrix_2d, 2)
        triangle = trianglesMatrix_2d{i,h};
        com = center_of_mass_matrix{i,h};
        for j = 1:3
            if (com(1) > d+(h-1)/2 && com(2) > L*sqrt(3)/2)
                trianglesMatrix_2d{i,h} = [];
                speciesMatrix(i,h) = 0;
            end
        end
    end
end

%% Calculate the final species vector
speciesVector_1 = zeros(num_triangles);

%Get everything in the first unit
for i = 1:2*c*(L+d)
    speciesVector_1(i) = speciesMatrix(mod(i-1, 2*c)+1, floor((i-1)/(2*c))+1);
end

speciesVector_1 = flip(speciesVector_1, 2);

% Get rid of the zero entries that correspond to empty triangles
speciesVector_1 = nonzeros(speciesVector_1)';

% Find the species of all the other triangles
speciesVector = speciesVector_1;
while length(speciesVector) < num_triangles
    speciesVector = [speciesVector, speciesVector_1(1:2*c*L), speciesVector_1(2*c*L+1:2*c*L+2*d^2)]; %it was for the second one: 2*c*L+1:2*L*(c+d)), but this fails for (3,2,1,1)
end


%% define the functions we use
function [xy] = hk_to_xy(h,k)
x_val = h + k/2;
y_val = k*sqrt(3)/2;
xy = [x_val,y_val];
end

function [hk] = xy_to_hk(x,y)
k_val = 2*y/sqrt(3);
h_val = x-k_val/2;
hk = [h_val,k_val];
end

function new_coords = unitPeriodic(ih, c,L,d)
new_coords = [mod(ih(1),c), mod(ih(2), L+d)];
end

function new_coords = unitPeriodic_xy(xy, c,L,d,x)
hk = xy_to_hk(xy(1), xy(2));
if hk(1) > c
    hk(1) = hk(1) - c;
end
if hk(1) < 0
    hk(1) = hk(1) + c;
end
if hk(2) > L+d %|| hk(2) < 0
    hk = hk - [x,0];
    hk(2) = hk(2) - (L+d);
end
if hk(2) < 0
    hk = hk + [x,0];
    hk(2) = hk(2) + L+d;
end
% periodic_hk = [mod(hk(1), c), mod(hk(2), L+d)]
new_coords = hk_to_xy(hk(1), hk(2));
end

function new_coords = unitCoords(triangle_coordinates, c,L,d)
periodic_location = [mod(triangle_coordinates(1),c), mod(triangle_coordinates(2), L+d)];
new_coords = [periodic_location(1) + periodic_location(2)/2, periodic_location(2)*sqrt(3)/2];
end

end
