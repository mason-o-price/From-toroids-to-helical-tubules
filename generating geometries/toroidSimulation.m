%% Simulation of a punctured tubule relaxing into a toroid or helical tubule.
% Stage 1: apply the force of the pseudo-bonds to close the holes and get near the final shape.
% Stage 2: merge the holes, and (if selected) maintain the cylindrical shape of the tubelets and push out the junctions.
% Stage 3: remove all pseudobonds and let the final structure relax.

%% ====================== Necessary Functions: ============================
% rot2d, available here: https://github.com/iswunistuttgart/matlab-tooling/blob/master/math/rot2d.m
% makeTubule
% removeHoles
% findBonds
% mergeHoles
% relaxEdges
% findUniqueTriangles
% plotEdges
% plotFaces
% makePeriodic
% findIndx
% findBindingAngles
% makeFigure
% findSpecies
% reorderTriangles
% reorderVertices
% simplifyTriangles
% simplifyConnections
% convertCoordinates
% map_tube_to_torus_func
% calculate_elastic_energy
% calculateStrain
% findInteractionMatrix
% measureOneAngle

%% ======================== Input Parameters ==============================
% clear all;

% Tubule parameters
T = 4;  % Lattice number (m,n) of the tube (m >= abs(n))

% Toroid parameters
numTubules = 4; % Number of tubule/junction sections around the toroid
tubuleHeight = 2; % Height of each tubule
junctionSideLength = 2; % Gives an LxL junction. This must be less than or equal to m/2
x = 1; % Relative position of the holes (between 0 and m-1). For a toroid, choose x = c-d-L/2
c = 1; % Offset for the hole locations (between 0 and m-1). This shouldn't change the initial structure (except for when the boundaries are merged). 
% NOTE: don't change c.

% Physical parameters
dt = 1e-1; % Time step per iteration
kConstTubule = 3; % Spring constant for the tubule section, this one should be high.
kConstJunction = 3; % Spring constant for the junction
kConstPseudoBond = 5e-1; % Spring constant for the pseudobonds, this one should be the weakest.
kConstJunction_pushout = 6e-2; % Spring constant for the force pushing out the junctions to keep them from buckling inward.
inflationStrength = 0;
kAngle = 0; % Bias strength of the angles
% Nice default values:
% dt = 1e-1;
% kConstTubule = 2;
% kConstJunction = 2;
% kConstPseudoBond = 2;
% kConstJunction_pushout = 6e-2;
% inflationStrength = 3e-2;
% kAngle = 2e-1;

% Simulation parameters
% (May need to tune for the simulation to work well)
iterationThreshold = 1.4e3; % Threshold for iteration
sampleFrequencyStage1 = 4; % Number of iterations between samples
sampleFrequencyStage2 = 12; % Number of iterations between samples
start_stage2 = 1e3; % Iteration when we close the holes
stop_forcing_tubelets = inf; % Iteration when we stop maintaining the tubes
stop_pushing_junctions = 1.2e3; % Iteration when we stop pushing out the junctions to keep them from buckling
start_biasing_angles = NaN; % Iteration when we start biasing the angles
stop_biasing_angles = NaN; % Iteration when we stop biasing the angles
angle_check_frequency = 10; % How often should we calculate the preferred angles
inflationStart = NaN; % Iteration to start pushing out on each vertex (i.e., inflating the mesh)
inflationEnd = NaN; % Iteration to stop inflating
pushout_junctions = 1; % Forcefully push out the junctions
maintain_tubes = 1; % Forcefully maintain a consistent diameter for the tubules
map_to_torus = 0; % 1 if toroid; 0 if helical tube. (Maps vertices onto the surface of a perfect torus).
twist_offset = 0; % Offset for the amount of twisting on the toroid
close_ends = 0; % 1 if toroid; 0 if helical tube.
% (Change when needed)
save_stretchingFigures = 0; % 1 if you want to save the figure of the colored edges based on stretching (png, svg, fig)
save_stretchingData = 0; % 1 if you want to save the stretching vector (mat, csv)
save_figures = 0; % 1 if you want to save the figures (png, svg, fig)
save_angles = 0; % 1 if you want to save the angles (mat, csv)
save_data = 1; % 1 if you want to save the coordinates, vertexConnections (mat, csv)
delete_empty_vertices = 0; % 1 if you need to use coordinates matrix; 0 if you only need the plots.
shouldAnimate = 1; % First animation; 1 = make a video; 0 = don't make a video
shouldAnimate_rotation = 1; % Second animation; 1 = make a video; 0 = don't make a video (for a rotating figure)
plot_energy = 0; % Make a plot of the elastic energy of the configuration as it evolves
% (Usually no need to change)
should_plot_markers = 1; % plot the line-segments indicating triangle orientation
addLight = 0; % Should we add light & shading
fixTop = 0; % Should we fix the top row of vertices in place
view_angle = [-29, 29]; %[-7, 8]; % [80, 20]; % What is the viewing angle? 3 = default, [Az, El] = anything
show_pseudoBonds = 1; % Plot the blue line-segments stitching the holes.
show_extraBonds = 0; % Plot the red line-segments keeping the tubules & junctions from buckling
plot_each_sample = 1; % Make a plot each time we take a sample
plot_stretching = 1; % 1 if you want to color the edge lengths based on the amount of stretching.
plot_angle_strain = 1; % 1 if you want to color the edges based on off-target angles w.r.t. the mean measured angle
shouldCalculateStrain = 1; % 1 if you want to calculate strain; 0 if not.

%% ============================ Preparation ===============================
% Find the total height of the tubule, (toroid's circumference)
totalHeight = numTubules*(tubuleHeight + junctionSideLength) + 1;

% Obtain the coordinates of the tubule vertices and the matrix of connections
[coordinates, vertexConnectivity, tube_radius] = makeTubule(m, 0, totalHeight);
% Format:
% coordinates matrix: row = vertex, columns 1,2,3,4,5 = i,h,x,y,z
% vertexConnectivity: cell matrix; row = i, column = h, cell = matrix of neighbors' (i,h) values.

% Make a matrix of the triangles
numTriangles_before = numTubules*m*(tubuleHeight + junctionSideLength)*2; %number of triangles before we cut the holes
numTriangles = numTubules*(2*tubuleHeight*m + 2*junctionSideLength^2); % Total number of triangles after we cut the holes
numSpecies = junctionSideLength^2 + m*tubuleHeight; %number of species
%trianglesMatrix = findUniqueTriangles(vertexConnectivity, coordinates, numTriangles_before, junctionSideLength, tubuleHeight, m, x); %Matrix of the triangles
% Format: trianglesMatrix = cell matrix; row n = nth triangle; only 1 column; cell = 3x2 matrix storing (i,h) values of the triangle vertices
trianglesMatrix = cell(numTriangles, 1);

% Find the triangles
counter = 1;
for k = 1:size(coordinates, 1) %go through all of the points
    point = coordinates(k,1:2);
    point = makePeriodic(point, m);
    i = point(1); %Find the i,h coordinates
    h = point(2);

    if ~isempty(vertexConnectivity{i,h}) %check if the point is deleted -> if the point's cell is not empty, continue...
        %Find the neighbors
        neighbors = vertexConnectivity{i,h};
        for j = 1:size(neighbors,1)
            n = neighbors(j, :); %This will give you the i,h values for this neighbor

            %If two vertices share another vertex in common --> they form a triangle
            j_neighbors = vertexConnectivity{n(1), n(2)};
            intersection = intersect(j_neighbors, neighbors, 'rows'); %Find the intersection of the cells.
            if ~isempty(intersection) %If the two vertices share a neighbor in common, draw a patch of a plane bounded by the 3 vertices
                for row = 1:size(intersection, 1)
                    final_vertex = intersection(row, :);
                    if ~isequal([h, h, h], [h, n(2), final_vertex(2)]) %Make sure that we aren't filling the cross-sections of a tubule, (only happens for a (3,0) tubule).
                        trianglesMatrix{counter, 1} = [i,h; n(1), n(2); final_vertex(1), final_vertex(2)];
                        counter = counter + 1;
                    end
                end
            end
        end
    end
end

% Plot the faces on the tubule
makeFigure(1); %set up the figure
view(view_angle); %set the view to be from the top-right
plotFaces(coordinates, trianglesMatrix, m, 0, 0, 0); %plot the faces.  color matrix: ones(numTriangles, 1), [1,1,1]);
if addLight
    light
end
shg %Display the figure

% Remove the holes
vertexConnectivity = removeHoles(m, x, c, numTubules, tubuleHeight, junctionSideLength, vertexConnectivity, coordinates, trianglesMatrix);

% Update the triangles matrix after we remove the holes.
trianglesMatrix = findUniqueTriangles(vertexConnectivity, coordinates, numTriangles, junctionSideLength, tubuleHeight, numTubules, m, x);

% Re-arrange the triangles in the triangles matrix
trianglesMatrix = reorderTriangles(trianglesMatrix,m,tubuleHeight,junctionSideLength,x,numTubules);

% Re-arrange the vertices within the triangles matrix
trianglesMatrix = reorderVertices(trianglesMatrix,m);

% Make a version of the triangles matrix that is not a cell matrix
trianglesIndx = simplifyTriangles(trianglesMatrix,m);

% Make a version of the vertex connections matrix that is not a cell matrix
%connectionsIndx = simplifyConnections(vertexConnectivity,coordinates,m);

% Mirror the vertices by reversing the sign of the z-coordinates
coordinates(:,5) = -coordinates(:,5);

% Find the right species associated with each triangle, and define a color matrix
speciesVector = findSpecies(m, tubuleHeight, junctionSideLength, x, numTriangles);
colorMatrix = videcolors(numSpecies); %  a different color for each species to color the faces. Thomas' color map :p
% colorMatrix = videcolors_for_422x;

% Set up the 2nd figure
makeFigure(2); 
view(view_angle); %set the view to be from the top-right

% Plot the faces
plotFaces(coordinates, trianglesMatrix, m, should_plot_markers, speciesVector, colorMatrix); 
if addLight
    light
end
shg % Display the figure

% Find the bonds across the holes, and extra bonds to make the tubules rigid
[Pcollection, Qcollection, pseudoBonds, extraBonds, junctionExtraBonds] = findBonds(numTubules, tubuleHeight, junctionSideLength, totalHeight, x, m, c);
% Format: 
% Pcollection: row = vertex, columns 1,2 = i,h
% Qcollection: row = vertex, columns 1,2 = i,h
% pseudoBonds: cell matrix; row = i, column = h, cell = vertex that the point (i,h) is pseudo-bound to.
% extraBonds: cell matrix; row = i, column = h, cell = vertices to which this point is connected.

% Find the default length of the extra bonds (see line 188)
defaultExtraBondLength = norm(coordinates(1, 3:5) - coordinates(1+floor(m/2), 3:5));
defaultLength2 = 2*defaultExtraBondLength;

% --------------- Map the tubule onto a torus ---------------
if map_to_torus
    coordinates = map_tube_to_torus_func(m, tubuleHeight, junctionSideLength, x, numTubules, tube_radius, coordinates, twist_offset);
    makeFigure(10); %set up the figure
    view(view_angle); %set the view to be from the top-right
    plotFaces(coordinates, trianglesMatrix, m, 1, speciesVector, colorMatrix);
    shg
end

% Generate an interaction matrix
interactionMatrix = findInteractionMatrix(m, tubuleHeight, junctionSideLength, x, numTriangles);
%% ======================= Run the simulation =============================
iteration = 1;

length = size(coordinates(:,1)); %Side-length of the square matrix below
kMatrix = zeros(length(1), length(1)); % Create a placeholder matrix of zeros to be filled with the k constants between connected vertices
% Format: Row = indx of vertex in coordinates matrix; Column = indx of neighbor's vertex

% Populate the k matrix
for k = 1:numTubules %Consider each unit
    %Start with the junction sections
    %Define the refernce corner a1
    a1 = [c+(k-1)*x, k*tubuleHeight + (k-1)*junctionSideLength]; %define the reference corner a1
    a1 = a1 + [0,1]; %offset the height to match MATLAB's indexing
    a1 = makePeriodic(a1, m);
    for u = 0:junctionSideLength
        for v = 1:junctionSideLength
            vertex = a1 + [u,v]; %Find a vertex inside of the junction
            vertex = makePeriodic(vertex, m);
            indx = (vertex(2)-1)*m+vertex(1); %Find the index
            neighbors = vertexConnectivity{vertex(1),vertex(2)}; %find the neighbors of this vertex
            for row = 1:size(neighbors,1)
                neighbor = neighbors(row, :);
                neighbor_indx = (neighbor(2)-1)*m+neighbor(1); %Find the index of the neighbors entry in the coordinates matrix
                kMatrix(indx, neighbor_indx) = kConstJunction; %Update the k matrix with the spring strength that we want for the connections within the junction
                kMatrix(neighbor_indx, indx) = kConstJunction; %Make sure the matrix is symmetric
            end
        end
    end

    %Now do the tubule sections
    for u = 1:tubuleHeight + 1 %Go through the horizontal rings of the tubule
        h = (k-1)*tubuleHeight + (k-1)*junctionSideLength + u; %Consider the height of this ring of vertices
        for i = 1:m %Cycle through all of the vertices in the ring
            indx = (h-1)*m + i;
            neighbors = vertexConnectivity{i,h}; %find the neighbors
            if ~isempty(neighbors)
                for row = 1:size(neighbors,1)
                    neighbor = neighbors(row, :);
                    neighbor_indx = (neighbor(2)-1)*m+neighbor(1); %Find the index of the neighbors entry in the coordinates matrix
                    kMatrix(indx, neighbor_indx) = kConstTubule; %Update the k matrix with the spring strength that we want for the tubule connections
                    kMatrix(neighbor_indx, indx) = kConstTubule; %Make sure the matrix is symmetric
                end
            end
        end
    end

    %Now do the pseudo-bonds
    for row = 1:size(Pcollection,1)
        p = Pcollection(row, :); %Find the points that are paired together
        q = Qcollection(row, :);
        p_indx = findIndx(p(1), p(2), m); %Find their indices in the coordinates matrix
        q_indx = findIndx(q(1), q(2), m);
        kMatrix(p_indx, q_indx) = kConstPseudoBond; %Update the k matrix
        kMatrix(q_indx, p_indx) = kConstPseudoBond; %Make sure it's symmetric
    end

    %Now do the extra bonds to make the tubules rigid
    for i = 1:size(extraBonds, 1)
        for h = 1:size(extraBonds, 2) %iterate through all of the entries in the extraBonds matrix
            if ~isempty(extraBonds{i,h}) %if the entry is not empty
                pnt = [i,h]; %find the point
                pnt = makePeriodic(pnt, m);
                pnt_indx = (pnt(2)-1)*m+pnt(1); %find the index of the point
                targets = extraBonds{pnt(1), pnt(2)}; %Find the points to connect to
                if ~isempty(targets) %make sure it's not empty
                    for row = 1:size(targets,1)
                        target = targets(row, :); %find a point to connect to
                        target_indx = (target(2)-1)*m+target(1);
                        kMatrix(pnt_indx, target_indx) = kConstTubule; %Update the kmatrix with a spring constant
                        kMatrix(target_indx, pnt_indx) = kConstTubule; %Do it for the transposed position as well.
                    end
                end
            end
        end
    end

    %Now do the junction extra bonds to make sure the junctions don't buckle
    for i = 1:size(junctionExtraBonds, 1)
        for h = 1:size(junctionExtraBonds, 2) %iterate through all of the entries in the extraBonds matrix
            if ~isempty(junctionExtraBonds{i,h}) %if the entry is not empty
                pnt = [i,h]; %find the point
                pnt = makePeriodic(pnt, m);
                pnt_indx = (pnt(2)-1)*m+pnt(1); %find the index of the point
                targets = junctionExtraBonds{pnt(1), pnt(2)}; %Find the points to connect to
                if ~isempty(targets) %make sure it's not empty
                    for row = 1:size(targets,1)
                        target = targets(row, :); %find a point to connect to
                        target_indx = (target(2)-1)*m+target(1);
                        kMatrix(pnt_indx, target_indx) = kConstJunction_pushout; %Update the kmatrix with a spring constant
                        kMatrix(target_indx, pnt_indx) = kConstJunction_pushout; %Do it for the transposed position as well.
                    end
                end
            end
        end
    end
end

%Prepare for making an animation by preallocating a movie structure.
mov(1:floor(start_stage2/sampleFrequencyStage1) + floor((iterationThreshold - start_stage2)/sampleFrequencyStage2)) = struct('cdata', [], 'colormap', []);
rotationMov(1:floor(start_stage2/sampleFrequencyStage1) + floor((iterationThreshold - start_stage2))) = struct('cdata', [], 'colormap', []);
frameCount = 0; 

stage = "Stage 1";
apply_pseudo_bonds = 1;

% Figure out where to start based on whether we are fixing the top row of points.
starting_pnt = 1;
if fixTop
    starting_pnt = m+1;
end

% make vectors storing the values of elastic energy and RMSD
el_energy = [];
RMSD_vector = [];
coordinates_record = coordinates; % make a copy of the previous version of the coordinates

% Add the forces
forceMatrix = zeros(1, 3); % Create a matrix for the forces. Format: row = vertex index, columns 1,2,3 = force components x,y,z
% Iterate the simulation
while (iteration < iterationThreshold)  
    %Go through all of the points. 
    for j = starting_pnt:size(coordinates,1)
        point = coordinates(j, :); % Find the point
        neighbors = vertexConnectivity{point(1), point(2)}; %find the neighbors of the point
        if ~isempty(neighbors) %if the point was not deleted, continue...
            indx = (point(2)-1)*m+point(1); %Find the index of the vertex in the coordinates matrix
            force = []; %Prepare a matrix storing the forces on each vertex
            %Apply force: neighbors
            for row = 1:size(neighbors, 1)
                neighbor = neighbors(row, :);
                neighbor_indx = (neighbor(2)-1)*m+neighbor(1);
                neighborXYZ = coordinates(neighbor_indx, 3:5);
                relativeVector = point(3:5) - neighborXYZ; %find the displacement between vertices
                displacement = relativeVector - relativeVector/norm(relativeVector); %Find the displacement of the relative vector from the rest-length
                force = [force; -1*kMatrix(indx, neighbor_indx)*displacement]; %Update the forces matrix
            end

            %Apply force: pseudo-bonds
            if apply_pseudo_bonds
                if ~isempty(pseudoBonds{point(1), point(2)}) %If this vertex has a pseudo-bond, add another force
                    q = pseudoBonds{point(1), point(2)};
                    q_indx = (q(2)-1)*m+q(1);
                    q_XYZ = coordinates(q_indx, 3:5);
                    displacement = point(3:5) - q_XYZ; % Find the displacement of the vertices from the rest length (which is 0 here)
                    force = [force; -1*kMatrix(indx, q_indx)*displacement]; % Update the forces matrix again with the pseudo-bond
                end
            end

            %Apply force: extra bonds
            if ( iteration < stop_forcing_tubelets ) && maintain_tubes
                if point(1) <= size(extraBonds, 1) && point(2) <= size(extraBonds, 2) % Make sure the point is not out of the bounds of the matrix
                    if ~isempty(extraBonds{point(1), point(2)}) % If this vertex has extra bonds, add another force
                        targets = extraBonds{point(1), point(2)}; % Find the set of vertics that this point is supposed to connect to
                        for row = 1:size(targets,1)
                            q = targets(row, :); % Choose a target point
                            q_indx = (q(2)-1)*m+q(1); % Find its index
                            q_XYZ = coordinates(q_indx, 3:5); % Find its xyz
                            relativeVector = point(3:5) - q_XYZ; % Find the vector connecting them
                            displacement = relativeVector - defaultExtraBondLength*relativeVector/norm(relativeVector); % Find the displacement of the vertices from the rest length (which is the diameter of the tubule here)
                            force = [force; -1*kMatrix(indx, q_indx)*displacement]; % Update the forces matrix again with the pseudo-bond
                        end
                    end
                end
            end

            %Apply force: junction extra bonds
            if ( iteration < stop_pushing_junctions ) && pushout_junctions
                if ( point(1) <= size(junctionExtraBonds, 1) ) && ( point(2) <= size(junctionExtraBonds, 2) ) % Make sure the point is not out of the bounds of the matrix
                    if ~isempty(junctionExtraBonds{point(1), point(2)}) % If this vertex has extra bonds, add another force
                        targets = junctionExtraBonds{point(1), point(2)}; % Find the set of vertics that this point is supposed to connect to
                        for row = 1:size(targets,1)
                            q = targets(row, :); % Choose a target point
                            q_indx = (q(2)-1)*m+q(1); % Find its index
                            q_XYZ = coordinates(q_indx, 3:5); % Find its xyz
                            relativeVector = point(3:5) - q_XYZ; % Find the vector connecting them
                            displacement = relativeVector - defaultLength2*relativeVector/norm(relativeVector); % Find the displacement of the vertices from the rest length (which is the diameter of the tubule here)
                            force = [force; -1*kMatrix(indx, q_indx)*displacement]; % Update the forces matrix again with the pseudo-bond
                        end
                    end
                end
            end
        end
        forceMatrix(j,:) = sum(force); % Add the total force on this point to the cell matrix of all the forces
    end

    inflateMatrix = zeros(size(forceMatrix));
    % Apply force: inflation
    if ( iteration < inflationEnd ) && ( iteration > inflationStart )
        for j = 1:size(trianglesIndx,1) % Go through all of the triangles. 
            triangle = trianglesIndx(j,:);
        
            % find the vertex indices
            v1Indx = triangle(1); % vertex 1
            v2Indx = triangle(2); % vertex 2
            v3Indx = triangle(3); % vertex 3
        
            v1XYZ = coordinates(v1Indx, 3:5); %find the xyz coordinates of v1
            v2XYZ = coordinates(v2Indx, 3:5); %find the xyz coordinates of v2
            v3XYZ = coordinates(v3Indx, 3:5); %find the xyz coordinates of v3
        
            % Determine the (counter-clockewise) oriented edge vectors on triangle A
            e1 = v1XYZ - v3XYZ;
            e2 = v2XYZ - v1XYZ;
            e3 = v3XYZ - v2XYZ;
        
            % Determine the unit face normal vector
            faceNorm = cross(e1, e2)/norm(cross(e1, e2));

            % Add the inflation force to each of the triangle's vertices
            inflateMatrix(v1Indx,:) = inflateMatrix(v1Indx,:) + faceNorm;
            inflateMatrix(v2Indx,:) = inflateMatrix(v2Indx,:) + faceNorm;
            inflateMatrix(v3Indx,:) = inflateMatrix(v3Indx,:) + faceNorm;
        end
    end
    
    if ( iteration < inflationEnd ) && ( iteration > inflationStart )
        % Normalize the vectors in the inflation matrix. Format: row = Fx, Fy, Fz.
        % Normalize: row --> (Fx, Fy, Fz)/||(Fx, Fy, Fz)||
        inflateMatrix = inflateMatrix.*1./(sum(inflateMatrix.*inflateMatrix,2).^0.5); % Take the square root to get the norms, and take the reciprocal so that we normalize the force of inflation
        forceMatrix = forceMatrix + inflationStrength*inflateMatrix; % Update the force matrix to include the inflation and weight the force vector by the inflation strength 
    end

    % Find average angles
    if (iteration >= start_biasing_angles) && ((mod(iteration, angle_check_frequency)==0) || iteration == start_biasing_angles)
        % Find the angles
        [~, angles] = findBindingAngles(coordinates, trianglesMatrix, speciesVector, numSpecies, m);
    end

    % Apply force: angle bias
    if ( iteration > start_biasing_angles ) && ( iteration <= stop_biasing_angles)
        % Go through each triangle
        for triangle_indx = 1:size(trianglesMatrix, 1)
            triangle = trianglesMatrix{triangle_indx, 1}; % find the vertices that make up this triangle. 
            species = speciesVector(triangle_indx); % find the species using the species vector
            for side = 1:3 % go through each side of this triangle
                neighborSpecies = floor((find(interactionMatrix(3*(species-1) + side, :) == 1) - 1)/3) + 1; % Find the species of the neighboring triangle
                neighborEdge = find(interactionMatrix(3*(species-1) + side, :) == 1) - 3*(neighborSpecies-1); % Find the side number of the neighbor
                possibleNeighborIndxs = find(speciesVector == neighborSpecies); % Find the possible indices of the neighbor
                
                hasNeighbor = 0; % Assume by default there is no neigbor on this edge
                % Find the exact index of the neighbor in the trianglesMatrix
                for j = 1:size(possibleNeighborIndxs, 2)
                    possibleNeighbor = trianglesMatrix{possibleNeighborIndxs(j), 1};
                    if size(intersect(triangle, possibleNeighbor, 'rows'), 1) == 2
                        neighbor = possibleNeighbor;
                        hasNeighbor = 1;
                        break
                    end
                end
                if hasNeighbor ~= 0
                    targetAngle = angles(species, neighborSpecies); % Find the target angle corresponding to this interaction
                    [measuredAngle, surfaceNormA, surfaceNormB] = measureOneAngle(triangle, neighbor, coordinates, m); % Measure the actual angle for this edge
                    dTheta = (measuredAngle - targetAngle)*pi/180; % find the difference in measured vs target angle, and convert to radians              
                    % Now apply this force to the vertex opposing this side. The opposing vertex to edge matching is:
                    % side 1 <--> vertex 2
                    % side 2 <--> vertex 3
                    % side 3 <--> vertex 1              
                    %            _                        _vertex 1         
                    %           /|\                      /|\            
                    %   side 2 / | \ side 1   <---->    / | \           
                    %         /_____\                  /_____\            
                    %          side 3          vertex 2       vertex 3    
                    vertexToChange = triangle(mod(side, 3) + 1, :); % find the coordinates of the vertex for the given side
                    neighborVertex = neighbor(mod(neighborEdge, 3) + 1, :); % find the neighbor's vertex to update
                    vertexToChangeIndx = findIndx(vertexToChange(1), vertexToChange(2), m); % find its index using the value from the triangle matrix
                    neighborVertexIndx = findIndx(neighborVertex(1), neighborVertex(2), m); % find the index of the neighbor's vertex.               
                    leverArm = sqrt(3)/2; % Assume a lever arm of sqrt(3)/2 for a perfect equilateral triangle. 
                    % Define the force so that it's pointing parallel to the surface norm of the triangle, and proportional to dTheta up to the bias strength kAngle
                    angleForceA = kAngle*dTheta^2*surfaceNormA/leverArm; 
                    angleForceB = kAngle*dTheta^2*surfaceNormB/leverArm; 
                    % sprintf("Force strength: %f", norm(angleForceA))
                    % disp("Final force:")
                    % disp(angleForce)  
                    % Update the force matrix at the corresponding vertex's location
                    forceMatrix(vertexToChangeIndx, :) = forceMatrix(vertexToChangeIndx, :) + angleForceA;
                    forceMatrix(neighborVertexIndx, :) = forceMatrix(neighborVertexIndx, :) + angleForceB;
                end
            end
        end
    end

    % Update all of the coordinates based on the forces
    for j = 1:size(coordinates, 1) %Go through all of the points
        point = coordinates(j,:);
        %if ~isempty(vertexConnectivity{point(1), point(2)}) %make sure this point isn't deleted
        force = forceMatrix(j,:);
        coordinates(j,3:5) = coordinates(j,3:5) + force*dt;
        %end
    end

    
    if ( mod(iteration, sampleFrequencyStage1)==0 && iteration<=start_stage2 ) || ( mod(iteration, sampleFrequencyStage2)==0 && iteration>start_stage2 )  %Take a sample at some frequency
        if iteration <= start_stage2
            fprintf(stage + ": [%.2f %%]\n", iteration/start_stage2*100) %Print the progress level
        else
            fprintf(stage + ": [%.2f %%]\n", (iteration-start_stage2)/(iterationThreshold-start_stage2)*100) %Print the progress level
        end

        if plot_each_sample
            clf(figure(3))
            figure(3);
            if plot_energy
                subplot(2,2,[1;3]); % make a figure with multiple plots. We use the left 2 for the 3d model.
                set(gcf, 'WindowStyle', 'normal')
                set(gcf, 'WindowState', 'maximized')
                
            end

            hold on; %make sure the plots don't overwrite each other
            set(gca,'visible','off'); %turn off the grid
            set(gcf,'Color','white'); %make the background white
            set(gcf, 'renderer','painters'); %make sure we can output vector file figures
            axis equal %set the axes to be equal size
            rotate3d(figure(3),'on');
            view(view_angle); %set the view to be from the top-right
            plotFaces(coordinates, trianglesMatrix, m, should_plot_markers, speciesVector, colorMatrix)
            if addLight
                light
            end
            %Plot the pseudo-bonds
            if apply_pseudo_bonds && show_pseudoBonds
                for i = 1:size(Pcollection,1)
                    startPnt = Pcollection(i,:);
                    startIndx = findIndx(startPnt(1), startPnt(2), m);
                    startXYZ = coordinates(startIndx, 3:5);
                    endPnt = Qcollection(i,:);
                    endIndx = findIndx(endPnt(1), endPnt(2), m);
                    endXYZ = coordinates(endIndx, 3:5);
                    lines = plot3([startXYZ(:,1)'; endXYZ(:,1)'],...
                        [startXYZ(:,2)'; endXYZ(:,2)'],...
                        [startXYZ(:,3)'; endXYZ(:,3)']);
                    set(lines, 'color', 'k', 'Linestyle', ':', 'LineWidth', 1); % [0,0.75,1]
                end
            end

            %Plot the extra bonds
            if show_extraBonds
                plotEdges(coordinates, extraBonds, m, 'r')
                plotEdges(coordinates, junctionExtraBonds, m, 'g')
            end

            % calculate and plot the elastic energy
            if plot_energy
                subplot(2,2,4)
                U = calculate_elastic_energy(vertexConnectivity, coordinates,m);
                el_energy = [el_energy, U];
                plot(1:iteration/size(el_energy,2):iteration, el_energy)
                title("Elastic energy")
                xlim([0,iteration])
                xlabel("Time (iterations)")
                ylabel("Potential energy")

                subplot(2,2,2)
                RMSD_vector = [RMSD_vector, calculate_RMSD(coordinates_record(:,3:5), coordinates(:,3:5))];
                coordinates_record = coordinates;
                plot(1:iteration/size(RMSD_vector,2):iteration, RMSD_vector)
                title("Root mean square deviation")
                xlim([0,iteration])
                xlabel("Time (iterations)")
                ylabel("RMSD")
            end

            %Add the frame to the animation
            frameCount = frameCount + 1;
            mov(frameCount) = getframe(gcf);
            drawnow;

            hold off;
        end
    end

    % Once we reach a certain point, start stage 2, and merge the holes
    if iteration == start_stage2
        apply_pseudo_bonds = 0; %remove the pseudo bonds
        maintain_tubules = 0; % stop pushing out tubule sections
        show_pseudoBonds = 0;

        % merge the holes
        [coordinates, vertexConnectivity, trianglesMatrix, kMatrix] = mergeHoles(coordinates, vertexConnectivity, trianglesMatrix, Pcollection, Qcollection, m, kMatrix, kConstTubule);

        if close_ends
            %%
            % Merge the open ends on the top and bottom
            top_end = []; % make a matrix storing the top & bottom vertices to be merged
            bottom_end = [];
            h = numTubules*(junctionSideLength + tubuleHeight)+1; %height along the top of the tiling
            %Top:
            for j = 1:m
                top_end = [top_end; [mod(j-x-1, m)+1, 1]];
            end
            %Bottom: 
            for j = 0:junctionSideLength-1
                bottom_end = [bottom_end; [mod((numTubules-1)*x+j,m)+1, h]];
            end
            for j = 1:m-2*junctionSideLength
                bottom_end = [bottom_end; [mod((numTubules-1)*x+2*junctionSideLength+j-1, m)+1, h-junctionSideLength]];
            end
            for j = 0:junctionSideLength-1
                bottom_end = [bottom_end; [mod((numTubules-1)*x,m)+1, h-junctionSideLength+j]];
            end
            [coordinates, vertexConnectivity, trianglesMatrix, kMatrix] = mergeHoles(coordinates, vertexConnectivity, trianglesMatrix, bottom_end, top_end, m, kMatrix, kConstTubule);
        end
        stage = "Stage 2"; %rename the stage
    end
    iteration = iteration + 1; % Increment the iteration
end

%% Calculate the binding angles
[bindingAngles_cell, bindingAngles] = findBindingAngles(coordinates, trianglesMatrix, speciesVector, numSpecies, m);
bindingAngles = round(bindingAngles, 4);

% save the angles
if save_angles
    save("DunlapToroids\bindingAngles\"+"("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_bindingAngles.mat", "bindingAngles")
    %csvwrite("DunlapToroids\bindingAngles\"+"("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_bindingAngles.csv",bindingAngles)
    dlmwrite("DunlapToroids\bindingAngles\"+"("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_bindingAngles.csv", bindingAngles, 'delimiter', ',', 'precision', 7); %using this function instead of csvwrite preserves more precision
end
%% Create an animation
if shouldAnimate == 1
    writerObj = VideoWriter("figures\toroid_("+m+","+tubuleHeight+","+junctionSideLength+","+x+")", 'MPEG-4');
    writerObj.FrameRate = 100; % set the images per second
    % writerObj.Quality = 100;
    open(writerObj); % open the video writer
    % write the frames to the video
    for i=1:size(mov, 2) - 1
        % convert the image to a frame
        frame = mov(i);
        %if isequal(size(frame.cdata,1:2), [writerObj.Width, writerObj.Height])
        if ~isempty(frame.cdata)
            writeVideo(writerObj, frame);
        end
        %end
    end
    % close the writer object
    close(writerObj);
end

%% Plot the final toroid
clf(figure(4));
makeFigure(4);
view(view_angle); %set the view 
if addLight
    light
end
%plotEdges(coordinates, vertexConnectivity, m, 'k');
plotFaces(coordinates, trianglesMatrix, m, should_plot_markers, speciesVector, colorMatrix);
%% Make a rotation animation
if shouldAnimate_rotation == 1
    for k = 1:200 %If you change this number, you also need to change the value in the rotationMov declaration.
        figure(4);
        hold on
        axis vis3d
        view( view_angle + [k*5, 0] );
        %Plot the pseudo-bonds
        if show_pseudoBonds
            for i = 1:size(Pcollection,1)
                startPnt = Pcollection(i,:);
                startIndx = findIndx(startPnt(1), startPnt(2), m);
                startXYZ = coordinates(startIndx, 3:5);
                endPnt = Qcollection(i,:);
                endIndx = findIndx(endPnt(1), endPnt(2), m);
                endXYZ = coordinates(endIndx, 3:5);
                lines = plot3([startXYZ(:,1)'; endXYZ(:,1)'],...
                    [startXYZ(:,2)'; endXYZ(:,2)'],...
                    [startXYZ(:,3)'; endXYZ(:,3)']);
                set(lines, 'color', 'b', 'Linestyle', '--');
            end
        end
        rotationMov(k) = getframe(gcf);
        drawnow;
    end
    writerObj1 = VideoWriter("figures\rotating_toroid_("+m+","+tubuleHeight+","+junctionSideLength+","+x+")", 'MPEG-4'); 
    writerObj1.FrameRate = 15; % set the images per second
    % writerObj1.Quality = 100;
    % open the video writer
    open(writerObj1);
    % write the frames to the video
    length = size(rotationMov);
    for i=1:length(2)-1
        % convert the image to a frame
        frame = rotationMov(i) ; 
        if ~isempty(frame.cdata) %,[writerObj1.Width, writerObj1.Height])
            writeVideo(writerObj1, frame);
        end
    end
    % close the writer object
    close(writerObj1);
end
%% Color the edges based on amount of stretching
if plot_stretching
    clf(figure(5));
    makeFigure(5);
    view(view_angle);
    stretching_vector = [];
    edges = []; %format: x1 x2,y1 y2,z1 z2 stretching
    for k = 1:size(coordinates,1) %go through all of the points
        point = coordinates(k,1:2);
        point = makePeriodic(point, m);
        i = point(1); %Find the i,h coordinates
        h = point(2);
        
        if h <= size(vertexConnectivity,2) && i <= size(vertexConnectivity, 1)
            if ~isempty(vertexConnectivity{i,h}) %check if the point is deleted -> if the point's cell is not empty, continue...
                %Find the neighbors
                neighbors = vertexConnectivity{i,h};
                for j = 1:size(neighbors,1)
                    n = neighbors(j, :); %This will give you the i,h values for this neighbor
                    nIndx = findIndx(n(1), n(2), m); %find the index in the coordinates matrix of xyz coordinates that corresonds to this point
                    neighborXYZ = coordinates(nIndx, 3:5); % neighbor's xyz coordinates
                    startXYZ = coordinates(k, 3:5); % starting point's xyz coordinates
                    stretch = norm(startXYZ - neighborXYZ) - 1; % find by how much the edges are stretched much the edges
                    stretching_vector = [stretching_vector, stretch];
                    edges = [edges; startXYZ(:,1), neighborXYZ(:,1), startXYZ(:,2), neighborXYZ(:,2), startXYZ(:,3), neighborXYZ(:,3)];
                end
            end
        end
    end

    %Plot the line segment
    max_stretch = max(abs(stretching_vector)); %find the maximum amount by which the edges are stretched
    min_stretch = -max_stretch; %min(stretching_vector(stretching_vector<0)); %find the maximum amount by which the edges are compressed
    % max_stretch = 0.0183;
    % min_stretch = -max_stretch; % -0.0223393;
    numColors = 200;
    cMap = colorcet('D04', N=numColors+1);
    
    for i = 1:size(edges,1)
        stretch = stretching_vector(i);
        %find the color to use
        colorIndx = int32(numColors+1-(max_stretch-stretch)/(max_stretch-min_stretch)*numColors);
        edgeColor = cMap(colorIndx,:);
        %plot the lines
         %{
            Format
            x1'; x2'
            y1'; y2'
            z1'; z2'
         %}
        lines = line([edges(i,1)'; edges(i,2)'],...
            [edges(i,3)'; edges(i,4)'],...
            [edges(i,5)'; edges(i,6)']);
        set(lines, 'color', edgeColor, 'LineWidth', 2)
   
    end
    
    %plot the faces (white, 0.5 FaceAlpha)
    for row = 1:size(trianglesMatrix, 1) 
        triangle = trianglesMatrix{row, 1}; % each row corresponds to a triangle
        XYZ = [];
        for j = 1:3
            XYZ = [XYZ; coordinates(findIndx(triangle(j,1), triangle(j,2), m), 3:5)];
        end
        patch(XYZ(:,1), XYZ(:,2), XYZ(:,3), [1,1,1], 'LineStyle', 'none', 'FaceAlpha', 0.75);%, colorMatrix(species,:));
    end
    
    %plot a color bar
    c=colorbar;
    tickCount = min_stretch : (max_stretch-min_stretch)/10 : max_stretch;
    c.TickLabels = num2cell(tickCount);
    colormap(cMap)
    title('edge length stretching')
    
    clf(figure(6))
    figure(6)
    histogram(stretching_vector);
    title('$\frac{\ell_{ij}}{\ell^{(0)}} - 1$', 'Interpreter','latex', 'FontSize', 20)
end

%% Plot the strain due to off-target angles
if plot_angle_strain
    clf(figure(7));
    makeFigure(7);
    view(view_angle);

    % prepare a colormatrix
    numColors = 200;
    cMapAngles = colorcet('D02', N=numColors+1);

    maxAngleDiffMatrix = [];
    minAngleDiffMatrix = [];
    for i = 1:size(trianglesMatrix, 1)
        triangle_A = trianglesMatrix{i,1}; % find the first triangle
        for j = 1:size(trianglesMatrix, 1)
            if i ~= j
                triangle_B = trianglesMatrix{j,1}; % find the second triangle
                if  size(intersect(triangle_A,triangle_B, 'rows'), 1) == 2
                    row = speciesVector(i);
                    col = speciesVector(j);
    
                    maxAngleDiffMatrix(row,col) = max(abs(bindingAngles_cell{row, col} - bindingAngles(row, col)));
                end
            end
        end
    end

    % Define the maximum off-target angle. NOTE: we take the minimum to be
    % the negative version so that the middle of the colorbar = 0 (i.e. so
    % that the values are symmetric). 
    maxAngleDiff = max(max(maxAngleDiffMatrix));
    minAngleDiff = -maxAngleDiff; %min(min(minAngleDiffMatrix));
    % maxAngleDiff = 58;
    % minAngleDiff = -58;

    for i = 1:size(trianglesMatrix, 1)
        triangle_A = trianglesMatrix{i,1}; % find the first triangle
        for j = 1:size(trianglesMatrix, 1)
            if i ~= j
                triangle_B = trianglesMatrix{j,1}; % find the second triangle
                if  size(intersect(triangle_A,triangle_B, 'rows'), 1) == 2 % if they share 2 vertices in common --> they are neighbors
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
                    % e3_A = v3_A_XYZ - v2_A_XYZ;
    
                    % Determine the (counter-clockwise) oriented edge vectors on triangle B
                    e1_B = v1_B_XYZ - v3_B_XYZ;
                    e2_B = v2_B_XYZ - v1_B_XYZ;
                    % e3_B = v3_B_XYZ - v2_B_XYZ;
    
                    % Determine the unit face normal vectors:
                    N_A = cross(e1_A, e2_A)/norm(cross(e1_A, e2_A));
                    N_B = cross(e1_B, e2_B)/norm(cross(e1_B, e2_B));
    
                    center_A = mean([v1_A_XYZ; v2_A_XYZ; v3_A_XYZ]); % Find the center of triangle A
                    center_B = mean([v1_B_XYZ; v2_B_XYZ; v3_B_XYZ]); % Find the center of triangle B
                    
                    d1 = norm(center_A - center_B); % find the initial distance between the triangle centers
                    
                    offset_A = center_A + 0.1*N_A; % follow the normal vector away from the center of triangle A
                    offset_B = center_B + 0.1*N_B; % follow the normal vector away from the center of triangle B
    
                    d2 = norm(offset_A - offset_B); % find the distance between triangles when measured from offset vectors in the direction of the face norms
    
                    if d2 < d1 % if the angle is negative (concave), then the offset vectors should become closer.
                        sign = -1; % set the angle to be negative.
                    else % otherwise, assume the angle is positive.
                        sign = 1;
                    end
    
                    % Calcualte the binding angle using atan2
                    measuredAngle = 180/pi*sign*abs(atan2(norm(cross(N_A, N_B)), dot(N_A, N_B)));
    
                    row = speciesVector(i); % species of triangle 1 = row in angles matrix
                    col = speciesVector(j); % species of triangle 2 = column in angles matrix
    
                    targetAngle = bindingAngles(row, col); % Retrieve the target angle from the binding angle matrix
    
                    angleDifference = measuredAngle - targetAngle; % Find the difference in the target and measured angles

                    % find the color that we should draw the edge
                    colorIndx = int32(numColors+1-(maxAngleDiff-angleDifference)/(maxAngleDiff-minAngleDiff)*numColors);
                    edgeColor = cMapAngles(colorIndx,:);
    
                    commonEdge = intersect(triangle_A, triangle_B, 'rows');
                    commonVertex1 = commonEdge(1,:);
                    commonVertex2 = commonEdge(2,:);
                    commonVertex1XYZ = coordinates(findIndx(commonVertex1(1), commonVertex1(2), m), 3:5);
                    commonVertex2XYZ = coordinates(findIndx(commonVertex2(1), commonVertex2(2), m), 3:5);
                    XYZ = [commonVertex1XYZ; commonVertex2XYZ];
                    
                    lines = line(XYZ(:,1), XYZ(:,2), XYZ(:,3)); % define the line for this edge
                    set(lines, 'color', edgeColor, 'LineWidth', 2) % draw the edge
                end
            end
        end
    end
    %plot the faces (white, 0.5 FaceAlpha)
    for row = 1:size(trianglesMatrix, 1) 
        triangle = trianglesMatrix{row, 1}; % each row corresponds to a triangle
        XYZ = [];
        for j = 1:3
            XYZ = [XYZ; coordinates(findIndx(triangle(j,1), triangle(j,2), m), 3:5)];
        end
        patch(XYZ(:,1), XYZ(:,2), XYZ(:,3), [1,1,1], 'LineStyle', 'none', 'FaceAlpha', 0.75);%, colorMatrix(species,:));
    end
    
    %plot a color bar
    c=colorbar;
    tickCount = minAngleDiff : (maxAngleDiff-minAngleDiff)/10 : maxAngleDiff;
    c.TickLabels = num2cell(tickCount);
    colormap(cMapAngles)
    title('average angle - measured angle')
end

%% Plot with the edges colored corresponding to the angles
if plot_angle_strain
    clf(figure(8));
    makeFigure(8);
    view(view_angle);

    % prepare a colormatrix
    numColors = 200;
    % cMapAngles = colorcet('D02', N=numColors+1);
    cMapAngles = videcolors(numColors+1);

    % Define the maximum off-target angle. NOTE: we take the minimum to be
    % the negative version so that the middle of the colorbar = 0 (i.e. so
    % that the values are symmetric). 
    maxAngle = 110;
    minAngle = -maxAngle; 

    for i = 1:size(trianglesMatrix, 1)
        triangle_A = trianglesMatrix{i,1}; % find the first triangle
        for j = 1:size(trianglesMatrix, 1)
            if i ~= j
                triangle_B = trianglesMatrix{j,1}; % find the second triangle
                if  size(intersect(triangle_A,triangle_B, 'rows'), 1) == 2 % if they share 2 vertices in common --> they are neighbors
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
                    % e3_A = v3_A_XYZ - v2_A_XYZ;
    
                    % Determine the (counter-clockwise) oriented edge vectors on triangle B
                    e1_B = v1_B_XYZ - v3_B_XYZ;
                    e2_B = v2_B_XYZ - v1_B_XYZ;
                    % e3_B = v3_B_XYZ - v2_B_XYZ;
    
                    % Determine the unit face normal vectors:
                    N_A = cross(e1_A, e2_A)/norm(cross(e1_A, e2_A));
                    N_B = cross(e1_B, e2_B)/norm(cross(e1_B, e2_B));
    
                    center_A = mean([v1_A_XYZ; v2_A_XYZ; v3_A_XYZ]); % Find the center of triangle A
                    center_B = mean([v1_B_XYZ; v2_B_XYZ; v3_B_XYZ]); % Find the center of triangle B
                    
                    d1 = norm(center_A - center_B); % find the initial distance between the triangle centers
                    
                    offset_A = center_A + 0.1*N_A; % follow the normal vector away from the center of triangle A
                    offset_B = center_B + 0.1*N_B; % follow the normal vector away from the center of triangle B
    
                    d2 = norm(offset_A - offset_B); % find the distance between triangles when measured from offset vectors in the direction of the face norms
    
                    if d2 < d1 % if the angle is negative (concave), then the offset vectors should become closer.
                        sign = -1; % set the angle to be negative.
                    else % otherwise, assume the angle is positive.
                        sign = 1;
                    end
    
                    % Calcualte the binding angle using atan2
                    measuredAngle = 180/pi*sign*abs(atan2(norm(cross(N_A, N_B)), dot(N_A, N_B)));
    

                    % find the color that we should draw the edge
                    colorIndx = int32(numColors+1-(maxAngle-measuredAngle)/(maxAngle-minAngle)*numColors); % Plot the color who's index is in proportion to the angle out of the total range of angles
                    edgeColor = cMapAngles(colorIndx,:);
    
                    commonEdge = intersect(triangle_A, triangle_B, 'rows');
                    commonVertex1 = commonEdge(1,:);
                    commonVertex2 = commonEdge(2,:);
                    commonVertex1XYZ = coordinates(findIndx(commonVertex1(1), commonVertex1(2), m), 3:5);
                    commonVertex2XYZ = coordinates(findIndx(commonVertex2(1), commonVertex2(2), m), 3:5);
                    XYZ = [commonVertex1XYZ; commonVertex2XYZ];
                    
                    lines = line(XYZ(:,1), XYZ(:,2), XYZ(:,3)); % define the line for this edge
                    set(lines, 'color', edgeColor, 'LineWidth', 2) % draw the edge
                end
            end
        end
    end
    %plot the faces (white, 0.5 FaceAlpha)
    for row = 1:size(trianglesMatrix, 1) 
        triangle = trianglesMatrix{row, 1}; % each row corresponds to a triangle
        species = speciesVector(row);
        XYZ = [];
        for j = 1:3
            XYZ = [XYZ; coordinates(findIndx(triangle(j,1), triangle(j,2), m), 3:5)];
        end
        patch(XYZ(:,1), XYZ(:,2), XYZ(:,3), [1,1,1], 'LineStyle', 'none', 'FaceAlpha', 1);%, colorMatrix(species,:));
    end
    
    %plot a color bar
    c=colorbar;
    tickCount = minAngle : (maxAngle-minAngle)/10 : maxAngle;
    c.TickLabels = num2cell(tickCount);
    colormap(cMapAngles)
    title('angle values')
end

%% Calculate the mean quadratic strain
if shouldCalculateStrain 
    [strain, lengthStrain, angleStrain] = calculateStrain(coordinates, vertexConnectivity, bindingAngles_cell, m);
    disp("Strain: " + strain);
    disp("( length: " + lengthStrain + ", angles: " + angleStrain + " )");
end

%% Get rid of all deleted vertices
if delete_empty_vertices
    for row = 1:size(coordinates,1)
        i = coordinates(row,1);
        j = coordinates(row,2);
        if isempty(vertexConnectivity{i,j})
            coordinates(row,:) = NaN(1,5); % if the vertex is supposed to be deleted, convert the row to a list of NaNs
        end
    end
    coordinates = coordinates(all(~isnan(coordinates),2),:);
end

%% save figures
if save_figures
    % save the final 3d structure as .fig, .svg, .png
    saveas(figure(4), "DunlapToroids\3Dfigures\("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_finalStructure3D.fig")
    saveas(figure(4), "DunlapToroids\3Dfigures\("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_finalStructure3D.svg")
    saveas(figure(4), "DunlapToroids\3Dfigures\("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_finalStructure3D.png")
end
%% Save data
if save_data
    % save the coordinates matrix as .mat
    save("DunlapToroids\coordinates\("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_coordinates.mat", "coordinates")
    writematrix(coordinates(:,3:5),"DunlapToroids\coordinates\("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_coordinates.csv") 
    
    % save the vertexConnectivity matrix as .mat
    save("DunlapToroids\vertexConnections\("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_vertexConnections.mat", "vertexConnectivity")
    writecell(vertexConnectivity,"DunlapToroids\vertexConnections\("+m+","+tubuleHeight+","+junctionSideLength+","+x+")_vertexConnections.csv") 
end
%% save the stretching figures & files
if save_stretchingFigures
    % save the stretch-based colored egdges as .fig, .svg, .png
    saveas(figure(5), "stretchingAnalysis\figures3D\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_coloredEdges.fig")
    saveas(figure(5), "stretchingAnalysis\figures3D\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_coloredEdges.svg")
    saveas(figure(5), "stretchingAnalysis\figures3D\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_coloredEdges.png")
    
    % save the histogram figure as .fig, .svg, .png
    saveas(figure(6), "stretchingAnalysis\histograms\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_stretchingDistribution.fig")
    saveas(figure(6), "stretchingAnalysis\histograms\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_stretchingDistribution.svg")
    saveas(figure(6), "stretchingAnalysis\histograms\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_stretchingDistribution.png")
    
end
if save_stretchingData
    % save the stretching vector as a .mat and .csv
    save("stretchingAnalysis\stretchingVectors\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_stretchingVector.mat", "stretching_vector")
    csvwrite("stretchingAnalysis\stretchingVectors\("+m+","+tubuleHeight+","+junctionSideLength+")_closed_stretchingVector.csv",stretching_vector)
end
%% Define a function to plot the neighborhood of an example point, (useful for trouble-shooting)
function [] = plotNeighborhood(i,h,m, coordinates, vertexConnectivity)
%i must be less than m
%h must be less than height
indx = findIndx(i, h, m);
plot3(coordinates(indx, 3), coordinates(indx, 4), coordinates(indx, 5), 'bo');

%find the neighbors
neighbors = vertexConnectivity{i,h};
for j = 1:size(neighbors(:,1))
    n = neighbors(j, :); %This will give you the i,h values for this neighbor
    n_indx = findIndx(n(1), n(2), m); %this will give you the index for the point in the total matrix
    plot3(coordinates(n_indx, 3), coordinates(n_indx, 4), coordinates(n_indx, 5), 'ro');
end

end

