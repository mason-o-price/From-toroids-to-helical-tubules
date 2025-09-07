function  plotEdges(coordinates, vertexConnections, m, edgeColor)
%PLOTEDGES Plot the set of edges for a curved tubule.
% 
% PLOTEDGES(coordinates, vertexConnections, m, edgeColor)
%
% INPUT:
%  coordinates # Matrix of vertex coordinates
%  vertexConnections # Set of vertex connections
%  m # Tubule circumference (T)
%  edgeColor # Color of the edges

view(3); %set the view to be from the top-right
set(gca,'visible','off'); %turn off the grid
set(gcf,'Color','white'); %make the background white
axis equal %set the axes to be equal size

for k = 1:size(coordinates(:,1)) %go through all of the points
    point = coordinates(k,1:2);
    point = makePeriodic(point, m);
    i = point(1); %Find the i,h coordinates
    h = point(2);
    
    if h <= size(vertexConnections,2) && i <= size(vertexConnections, 1)
        if ~isempty(vertexConnections{i,h}) %check if the point is deleted -> if the point's cell is not empty, continue...
            %Find the neighbors
            neighbors = vertexConnections{i,h};
            for j = 1:size(neighbors(:,1))
                n = neighbors(j, :); %This will give you the i,h values for this neighbor
                nIndx = findIndx(n(1), n(2), m); %find the index in the coordinates matrix of xyz coordinates that corresonds to this point
                neighborXYZ = coordinates(nIndx, 3:5); % neighbor's xyz coordinates
                startXYZ = coordinates(k, 3:5); % starting point's xyz coordinates
                
                %Plot the line segment
                lines = line([startXYZ(:,1)'; neighborXYZ(:,1)'],...
                    [startXYZ(:,2)'; neighborXYZ(:,2)'],...
                    [startXYZ(:,3)'; neighborXYZ(:,3)']);
                set(lines, 'color', edgeColor);
            end
        end
    end
end


