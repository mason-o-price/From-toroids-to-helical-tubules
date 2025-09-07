function matrix = findInteractionMatrix(T,L,D,R, numTriangles)
%FINDINTERACTIONMATRIX find the interaction matrix for given input parameters (T,L,D,R)
% 
% matrix = findInteractionMatrix(T,L,D,R, numTriangles)
% 
% INPUT:
%  T # Tubelet circumference
%  L # Tubelet length
%  D # Junction side-length
%  R # Rolling parameter
%  numTriangles # Number of triangle types
%
% OUTPUT:
%  matrix # Interaction matrix

%% Preparation
% Update variable names
c = T;
d = D;
x = R;
speciesVector = findSpecies(c, L, d, x, numTriangles);
num_species = c*L + d^2;

matrix = zeros(3*num_species, 3*num_species);
%% Find interactions
%interactions on side 1 & 2 on the interior of the tube
for i = 1:c*L %only need to go through the first half of the tubule segment
    if mod(i, 2*c) ~= 0 %make sure that we don't consider the triangles on the right edge
        species_1 = speciesVector(i);
        species_2 = speciesVector(i+1); %find the species of the triangle to the right
        if mod(i,2) ~= 0 %if odd, triangle is rightside up, binds on side 1 to the right
            interactionSide = 1;
        elseif mod(i,2) == 0 %if even, triangle is upside down, binds on side 2 to the right
            interactionSide = 2;
        end
        matrix(3*(species_1-1)+interactionSide, 3*(species_2-1)+interactionSide) = 1; %append the interaction to the matrix
        matrix(3*(species_2-1)+interactionSide, 3*(species_1-1)+interactionSide) = 1; %make it symmetric
    end
end
%interactions on side 2 on the seam of the tube
for k = 1:L
    interactionSide = 2;
    species_1 = speciesVector(2*c*(k-1)+1); %find the first species, this will bind along side 2
    species_2 = speciesVector(2*c*k); %find the second species, this will bind also along side 2
    matrix(3*(species_1-1)+interactionSide, 3*(species_2-1)+interactionSide) = 1; %append the interaction to the matrix
    matrix(3*(species_2-1)+interactionSide, 3*(species_1-1)+interactionSide) = 1; %make it symmetric
end

%interactions between the top of the tubule and the bottom of the junction
for k = 1:d
    species_1 = speciesVector(2*c*(L-1)+2*k);
    species_2 = speciesVector(2*c*L+2*k-1);
    interactionSide = 3;
    matrix(3*(species_1-1)+interactionSide, 3*(species_2-1)+interactionSide) = 1; %append the interaction to the matrix
    matrix(3*(species_2-1)+interactionSide, 3*(species_1-1)+interactionSide) = 1; %make it symmetric
end

%interactions on side 3 between consecutive layers in the tubule
if L > 1
    for k = 1:c*L %go through the first half of the tubule
        if mod(k,2) == 0 %take only the triangles pointing down
            species_1 = speciesVector(k);
            species_2 = speciesVector(k + 2*c-1);
            interactionSide = 3;
            matrix(3*(species_1-1)+interactionSide, 3*(species_2-1)+interactionSide) = 1; %append the interaction to the matrix
            matrix(3*(species_2-1)+interactionSide, 3*(species_1-1)+interactionSide) = 1; %make it symmetric
        end
    end
end

%interactions on sides 1 & 2 on the interior of the junction
for k = 1:d^2 %go through the first half of the triangles
    if mod(2*c*L+k, 2*c*L+2*d) ~= 0 %don't consider the triangle along the right edge of the junction
        species_1 = speciesVector(2*c*L+k);
        species_2 = speciesVector(2*c*L+k+1);
        if mod(k,2) ~= 0 %if odd, triangle is rightside up, binds on side 1 to the right
            interactionSide = 1;
        end
        if mod(k,2) == 0 %if even, triangle is upside down, binds on side 2 to the right
            interactionSide = 2;
        end
        matrix(3*(species_1-1)+interactionSide, 3*(species_2-1)+interactionSide) = 1; %append the interaction to the matrix
        matrix(3*(species_2-1)+interactionSide, 3*(species_1-1)+interactionSide) = 1; %make it symmetric 
    end
end

%interactions along the open sides of the junction (only need to do the right or left by symmetry. Here, do right so we don't have to account for x)
for i = 1:d
    species_1 = speciesVector(2*c*(L-1) + 2*d + 2*i); %along the top of the tubule
    interactionSide1 = 3;
    species_2 = speciesVector(2*c*L+2*d*i); %along the right of the junction
    interactionSide2 = 2;
    matrix(3*(species_1-1)+interactionSide1, 3*(species_2-1)+interactionSide2) = 1; %append the interaction to the matrix
    matrix(3*(species_2-1)+interactionSide2, 3*(species_1-1)+interactionSide1) = 1; %make it symmetric 
end

%now do interactions from the top to the bottom of a hole (i.e. between consecutive tubules)
for k = 1:c-2*d
    interactionSide = 3;
    species_1 = speciesVector(2*c*(L-1) + 4*d + 2*k); % first species is the species of the triangle at this location
    species_2 = speciesVector(2*c*L - 2*(k-1)); % second species ...
    matrix(3*(species_1-1)+interactionSide, 3*(species_2-1)+interactionSide) = 1; %append the interaction to the matrix
    matrix(3*(species_2-1)+interactionSide, 3*(species_1-1)+interactionSide) = 1; %Make sure the matrix is symmetric
end

%do side 3 in the interior of the junction
for i = 1:d^2
    if mod(i,2) == 0 %only do it for triangles pointing down (the interactions of the triangles pointing up should be taken care of by symmetry)
        species_1 = speciesVector(2*c*L + i);
        species_2 = speciesVector(2*c*L + i + 2*d-1);
        interactionSide = 3;
        matrix(3*(species_1-1)+interactionSide, 3*(species_2-1)+interactionSide) = 1; %append the interaction to the matrix
        matrix(3*(species_2-1)+interactionSide, 3*(species_1-1)+interactionSide) = 1; %make it symmetric 
   end
end
%concern: is the species vector meant for the 2d tiling or for the coordinates matrix?

% problems so far:
% - across the hole: side 3 along the tube -> side 2 along the junction

end