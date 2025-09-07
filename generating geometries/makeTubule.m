function [allCoordinates, vertexConnectivity, R] = makeTubule(m,n,height)
%MAKETUBULE Create coordinates and vertex connections matrices for a tube
% 
% [allCoordinates, vertexConnectivity, R] = MAKETUBULE(m,n,height)
% 
% INPUT:
%  m # first lattice number
%  n # second lattice number
%  height # desired height of the tubule
%
% OUTPUT:
%  allCoordinates # coordinates matrix. Format: [i,h,x,y,z] where i,h are the indices of the given vertex, and xyz are its coordinates. 
%  vertexConnectivity # cell matrix of vertex connections. Format: entry at {i,h} is the matrix of neighboring vertex indices. 
%  R # tubule radius
%% Solve for initial variables
%clf(figure(1))

height = height + 1; %We don't use all the vertices, so the actual height is less than what we would want.

%rotate the tube by this much from original
rotateby = 0; %degrees

absn = abs(n);
if n == 0 
    chiral = 1;
else
    chiral = n/abs(n); %(+1 for right handed, -1 for left handed)
end

%obtain seam (dipping) angle
%this is absolute angle of vector a1 to xy surface
phi = atan(sqrt(3) * absn / (2 * m + absn)); 

%a1, a2, a3 are the 3 vector component in 2D coordinate
%rotate a1 to get a2 and a3
a1 = [cos(phi); chiral * sin(phi)];
a2 = rot2d(-chiral*60) * a1;
a3 = rot2d(-chiral*120) * a1;

%look at Hayakawa et al. 2022 supplementary for reasoning behind the
%equations

R = 1/(2*sin(pi/m));
t1 = 2*pi/m;
x0 = [t1, t1/2, R]; %pick parameter starting point. A magically good value is often x0 = [0.7, 0, 1.8], or sometimes x0 = [1.3, 0, 4.5]. Unless you have a large radius tubule, in which case use a larger R. (last time: [1.6, 0, R])
%x0 = [0.9, 0, 6];
%Format: [theta1, theta2, R]
%theta1 is the angular component covered by one side along m
%theta2 is the angular component covered by one side along n
%Example: (3,0) = [1.6,0,0.5]

%options = optimoptions('fsolve', 'MaxFunctionEvaluations', 1e3, 'MaxIterations',1e3);

%x1 = fsolve(fun,x0, options);
x1 = x0; % in general we would have to solve the equations, but since we always do achiral tubules we can predict the parameters a priori

%if the values contain imaginary, use a different starting parameter
if abs(imag(x1(1))) > 10^-5
    error('Please choose better x0');
else
    theta1 = real(x1(1));
    theta2 = real(x1(2));
    theta3 = theta2-theta1;
    R = real(x1(3));
end

%% Unit vectors in cylindrical coordinate [r,theta,z]

%We take the unit vectors and project them onto a cylinder.
a1cyl = [0,theta1,a1(2)];
a2cyl = [0,theta2,a2(2)];
a3cyl = [0,theta3,a3(2)];

%Make a matrix of zeros as placeholders for the points.
numPoints = m*(height+1); %Set the number of points
plotpointCyl = zeros(numPoints,5);


%In all cases, we start the lattice at a right-most vertex, and add more to
%the left. Then, if the tubule is right-handed, that means we need to limit
%the bottom. If the tubule is left-handed, that means we need to limit the
%top. If the tubule is achiral, we limit the bottom (i.e., it grows
%downward).
if chiral == 1
    downlimit = height-1;
    uplimit = 0;
else
    downlimit = 0;
    uplimit = height -1;
end

for h = 0 : height
    for i = 0 : m-1

        %record m,n,R,theta,z of the point in the matrix
        plotpointCyl(h*m + i + 1,1) = i;
        plotpointCyl(h*m + i + 1,2) = h;
        plotpointCyl(h*m + i + 1,3) = R;

        %For each row in the matrix plotpointCyl:
        %The 1st entry is the position along a horizontal strip (between 0 & m+1)
        %The 2nd entry is the height of the horizontal strip (between 0 & height)
        %The 3rd entry is the radius. This should be the same for all.
        %The 4th entry is the angle of the other vertex in the triangle.
        %We only need two vertices, because the 3rd is captured by either
        %the previous or the next triangle
        %The 5th entry is the z position of the other vertex.

        %track theta and z and just keep on adding to it
        if i == 0
            theta = h * theta2;
            z = h * a2(2);
        else
            theta = theta + theta1;
            z = z + a1(2);
        end

        plotpointCyl(h*m + i + 1,4) = theta + rotateby*pi()/180;
        plotpointCyl(h*m + i + 1,5) = z;
    end
end

%% convert cylindrical to xyz coordinates

plotpointXYZ = zeros(numPoints,5); %Create a matrix of zeros the same size as the plotpointCyl matrix
plotpointXYZ(:,1:2) = plotpointCyl(:,1:2); %The first two columns are the lattice positions of each vertex, (i,h)
plotpointXYZ(:,3) = R * cos(plotpointCyl(:,4)); %The third column is the x-coordinates
plotpointXYZ(:,4) = R * sin(plotpointCyl(:,4)); %The fourth column is the y-coordinates
plotpointXYZ(:,5) = plotpointCyl(:,5); %The fifth column is the z-coordinates.

%% Find the vertices 
%Define vectors with the x,y,z coordinates
XYZ = plotpointXYZ(:, 3:5);

vLength = n+1:m*(height-1); %Restrict the vertices that we plot. This helps match them up with the line-segments below.
verticesToPlot = plotpointXYZ(vLength, 3:5);
allCoordinates = plotpointXYZ(vLength, :); %This is the first output of the function
allCoordinates(:, 1:2) = allCoordinates(:, 1:2) + 1; %Adjust the index in the first column:


%Plot the points
%plot3(verticesToPlot(:,1), verticesToPlot(:,2), verticesToPlot(:,3), 'k.'); 

%% Plot the line segments

% 1 -> 2 (black)
% 1 -> 3 (red)
% 2 -> 3 (blue)

%To prevent the top of the tubule from connecting to the bottom, we need to
%restrict the set of vertices we connect. 
%notation: length_ab = interval of vertices to connect between vertices a & b
length_12 = m+n+1:m*(height-1); %Set the length
length_23 = 1:m*(height-1); 

%Set the starting vertices for the line-segments
XYZstart_12 = XYZ(length_12, 1:3);
XYZstart_23 = XYZ(length_23, 1:3);

%Make a matrix with shifted entries to give a set of end-points to draw the
%lines between.
XYZend_12 = circshift(XYZ, m);
XYZend_13 = circshift(XYZ, m-1);
XYZend_23 = circshift(XYZ, 1);

%Limit the length of the vectors, to match the length of the start matrix
XYZend_12 = XYZend_12(length_12, 1:3);
XYZend_13 = XYZend_13(length_12, 1:3);
XYZend_23 = XYZend_23(length_23, 1:3);

%Reassign the seam edges
XYZend_13(2*m:m:end, :) = XYZend_13(m:m:end-m, :);
XYZend_23(1:m:end-m, :) = XYZend_23(m+1:m:end, :);
%We need to redefine the location that the final line-segments along each 
%ring get mapped to. The syntax (m:m:end, :) will redefine the entries
%starting at row m, skip m rows, and stop at the end of the matrix, 
%for all columns.

%Address the final missing edges
b = length(XYZend_23(:,1));
XYZend_23(b-m+1, :) = XYZ(b, :);
XYZend_13(m, :) = XYZ(1, :);

%% Define the Vertex Connectivity Matrix to find the neighbors of each point

vertexConnectivity = cell(m, height); %create a cell matrix
indx = plotpointXYZ(vLength, 1:2); %defie a matrix of all the indices (i,h)
indx = indx + 1;

startIndx_12 = indx(length_12, 1:2); %define the starting indices on edges 12 and 13
startIndx_23 = indx(length_23, 1:2); %do the same for edge 23

%Set the index for the end-points
endIndx_12 = circshift(indx, m);
endIndx_13 = circshift(indx, m-1);
endIndx_23 = circshift(indx, 1);

%limit the lengths as before
endIndx_12 = endIndx_12(length_12, 1:2);
endIndx_13 = endIndx_13(length_12, 1:2);
endIndx_23 = endIndx_23(length_23, 1:2);

%Reassign the seam edges
endIndx_13(2*m:m:end, :) = endIndx_13(m:m:end-m, :);
endIndx_23(1:m:end-m, :) = endIndx_23(m+1:m:end, :);

%Address the final missing edges
b = length(endIndx_23(:,1));
endIndx_23(b-m+1, :) = indx(b, :);
endIndx_13(m, :) = indx(1, :);

%Do side 12 & 13
for j = 1:size(startIndx_12, 1)
    %Define the index values for side 12
    i = startIndx_12(j, 1);
    h = startIndx_12(j, 2);
    vertexConnectivity{i,h} = [vertexConnectivity{i,h}; endIndx_12(j, 1:2)]; %add them to the matrix

    %Do the same for side 13
    i = startIndx_12(j, 1);
    h = startIndx_12(j, 2);
    vertexConnectivity{i,h} = [vertexConnectivity{i,h}; endIndx_13(j, 1:2)];
end

%Do side 23
for j = 1:size(startIndx_23, 1)
    i = startIndx_23(j, 1);
    h = startIndx_23(j, 2);
    vertexConnectivity{i,h} = [vertexConnectivity{i,h}; endIndx_23(j, 1:2)];
end

%Define a matrix that is exactly the vertexConnectivity matrix at this
%state. That is, a matrix of the unique connections between the vertices,
%so that we can draw the edges without double-counting.
%uniqueConnections = vertexConnectivity;

%% Fill out the rest of the vertexConnectivity matrix

%Now we need to do the other direction, i.e. instead of 12 we do 21 etc. We
%do this by reversing the start and end matrices that we use.
%Do side 21 & 31
for j = 1:size(endIndx_12, 1)
    %Define the index values for side 21
    i = endIndx_12(j, 1);
    h = endIndx_12(j, 2);
    vertexConnectivity{i,h} = [vertexConnectivity{i,h}; startIndx_12(j, 1:2)]; %add them to the matrix
    
    %Do the same for side 31
    i = endIndx_13(j, 1);
    h = endIndx_13(j, 2);
    vertexConnectivity{i,h} = [vertexConnectivity{i,h}; startIndx_12(j, 1:2)];
end

%Do side 32
for j = 1:size(endIndx_23, 1)
    i = endIndx_23(j, 1);
    h = endIndx_23(j, 2);
    vertexConnectivity{i,h} = [vertexConnectivity{i,h}; startIndx_23(j, 1:2)];
end

%% Relax the edges
%[allCoordinates, placeholder] = relaxEdges(allCoordinates, vertexConnectivity, 1e-1, m, 1, 1e3, 3e2, 0, 0, "Stage 1");

%% Plot the edges & the faces
%plotEdges(allCoordinates, vertexConnectivity, m, 'k')

%% Define a function to solve for theta1,2,3 and z
function F = solveTube(x,m,n)
%theta is the two value vectors

theta1 = x(1);
theta2 = x(2);
r = x(3);

theta3 = theta2-theta1;

z1 = sqrt(1-4*r^2*sin(theta1/2)^2);
z2 = sqrt(1-4*r^2*sin(theta2/2)^2);

z3 = z2-z1;

F(1) = m * z1 - n * z3;
F(2) = m * theta1 - n * theta3 - 2*pi;
F(3) = 4*r^2*sin((theta2-theta1)/2)^2 + (z2-z1)^2 -1;
end

end

