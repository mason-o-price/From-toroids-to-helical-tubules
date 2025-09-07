function [coordinates, mov] = relaxEdges(coordinates, vertexConnectivity, kConst, m, iteration, iterationThreshold, sampleFrequency, shouldAnimate, mov, stage)
%RELAXEDGES Minimize elastic energy between edges in a triangular mesh
%
% [coordinates, mov] = RELAXEDGES(coordinates, vertexConnectivity, kConst, m, iteration, iterationThreshold, sampleFrequency, shouldAnimate, mov, stage)
%
% INPUT:
%  coordinates # Matrix of veretx coordinates
%  vertexConnectivity # Matrix of vertex connections
%  kConst # Spring constant
%  m # Tubelet circumference (T)
%  iteration # Initial time step
%  iterationThreshold # Max number of time steps
%  sampleFrequency # Number of time steps between plots
%  shouldAnimate # Boolean whether we should plot the relaxation
%  mov # Movie structure
%  stage # Current stage number (1,2,3)
%
% OUPUT:
%  coordinates # updated coordinates
%  mov # updated movie structure

while (iteration < iterationThreshold)
    for j = 1:size(coordinates(:,1)) %Go through all of the points. 
        point = coordinates(j, :); % Find the point
        neighbors = vertexConnectivity{point(1), point(2)}; %find the neighbors of the point
        if ~isempty(neighbors) %if the point was not deleted, continue...
            indx = findIndx(point(1), point(2), m); %Find the index of the vertex in the coordinates matrix
            force = []; %Prepare a matrix storing the forces on each vertex
            %Force: neighbors
            for row = 1:size(neighbors(:,1))
                neighbor = neighbors(row, :);
                neighbor_indx = findIndx(neighbor(1), neighbor(2), m);
                neighborXYZ = coordinates(neighbor_indx, 3:5);
                relativeVector = point(3:5) - neighborXYZ; %find the displacement between vertices
                displacement = relativeVector - relativeVector/norm(relativeVector); %Find the displacement of the relative vector from the rest-length
                force = [force; -1*kConst*displacement]; %Update the forces matrix
            end
        end
        forceMatrix{point(1), point(2)} = sum(force); %Add the total force on this point to the cell matrix of all the forces
        % Update all of the coordinates based on the forces
    end
    for j = 1:size(coordinates(:,1)) %Go through all of the points
        point = coordinates(j, :);
        if ~isempty(vertexConnectivity{point(1), point(2)}) %make sure this point isn't deleted
            force = forceMatrix{point(1), point(2)};
            coordinates(j, 3:5) = coordinates(j, 3:5) + force;
        end
    end
    if mod(iteration, sampleFrequency) == 0 %Take a sample at some frequency
        fprintf(stage + ": [%.2f %%]\n", iteration/iterationThreshold*100) %Print the progress level

        if shouldAnimate == 1 
            clf(figure(3)); %clear the figure
            figure(3); %make a figure
            hold on %make sure it doesn't overwrite anything
            plotEdges(coordinates, vertexConnectivity, m, 'k'); %plot the edges

            %Add the frame to the animation
            mov(iteration/sampleFrequency) = getframe(gcf);
            drawnow;
        end
    end

    iteration = iteration + 1; % Increment the iteration
end

end

