function [Pcollection, Qcollection, pseudoBonds, extraBonds, junctionExtraBonds] = findBonds(numTubules, tubuleHeight, junctionSideLength, totalHeight, x, m, c)
%FINDBONDS Find the set of bonds between vertices
%
% [Pcollection, Qcollection, pseudoBonds, extraBonds, junctionExtraBonds] = findBonds(numTubules, tubuleHeight, junctionSideLength, totalHeight, x, m, c)
% 
% INPUT:
%  numTubules # number of tubelet segments (N)
%  tubuleHeight # height of the tubelet segments (L)
%  junctionSideLength # side-length of the junction regions (D)
%  totalHeight # total number of triangles along curved tubule axis N*(L+D)
%  x # shifting parameter between cuts in the holey tiling (R)
%  m # tubelet circumference (T)
%  c # lattice offset (assume 0)
%
% OUTPUT:
%  Pcollection # Matrix of starting vertices in stitching pairs.
%  Qcollection # Matrix of terminal vertices in stitching pairs.
%  pseudoBonds # Cell matrix summarizing the stitches across holes.
%  extraBonds # Additional bonds to keep the tubelet circular.
%  junctionExtraBonds # Pseudo-bonds to push out the junction. 

Pcollection = []; 
Qcollection = []; 
pseudoBonds = cell(m, totalHeight); % Cell matrix. Format: row = i, column = h, cell = vertex that the point (i,h) is pseudo-bound to.
extraBonds = cell(m, numTubules*tubuleHeight); % Cell matrix. Format: row = i, column = h, cell = vertex that the point (i,h) is pseudo-bound to.

% Fill up P and Q
for k = 1:numTubules
    %Define the corners
    a1 = [c+(k-1)*x, k*tubuleHeight + (k-1)*junctionSideLength]; %define the reference corner a1
    a1 = a1 + [0,1]; %offset the height to match MATLAB's indexing
    a1 = makePeriodic(a1, m);
    a2 = a1 + [junctionSideLength, 0]; %Define the first 5-fold corner
    a2 = makePeriodic(a2, m); %Implement the periodic boundary
    a4 = a1 + [0, junctionSideLength]; %Define the second 5-fold corner
    a4 = makePeriodic(a4, m);

    %Find all of the points from the junction that we are going to pair
    for j = 1:junctionSideLength
        %Start with the a2 corner
        p_a2 = a2 + [j, 0]; %find the first point
        p_a2 = makePeriodic(p_a2, m);
        q_a2 = a2 + [0, j]; %find the second point
        q_a2 = makePeriodic(q_a2, m);

        %Add them to the P/Qcollection matrices
        Pcollection = [Pcollection; p_a2]; %a2 corner first
        Qcollection = [Qcollection; q_a2];

        if a4(2) ~= totalHeight %make sure we don't stitch the very bottom corner to an empty vertex.
            %Now do the a4 corner
            p_a4 = a4 - [j, 0]; %find the first point 
            p_a4 = makePeriodic(p_a4, m);
            q_a4 = a4 - [0, j]; %find the second point
            q_a4 = makePeriodic(q_a4, m);      
            if q_a2 ~= p_a4 %make sure we're not overlapping in the middle
                Pcollection = [Pcollection; p_a4];
                Qcollection = [Qcollection; q_a4];
            end

        end
    end   
    %Find all of the points across the tubules that we are going to pair
    if junctionSideLength < m/2
        for i = 1:m-2*junctionSideLength
            if a4(2) ~= totalHeight %make sure we don't stitch the very bottom corner to an empty vertex.
                %add new points to stitch
                p = a2 + [i + junctionSideLength, 0];
                p = makePeriodic(p, m);
                q = p + [-junctionSideLength, junctionSideLength];
                q = makePeriodic(q, m);
                
                %Add them to the matrices
                Pcollection = [Pcollection; p];
                Qcollection = [Qcollection; q];
            end
        end
    end
end

% Populate the pseudoBonds matrix
for row = 1:size(Pcollection,1)
    p = Pcollection(row, :); % Find the first point
    q = Qcollection(row, :); % Find the second point which corresopnds to p
    pseudoBonds{p(1), p(2)} = [q(1), q(2)]; % Update the pseudoBonds matrix with the vertices that get paired.
    pseudoBonds{q(1), q(2)} = [p(1), p(2)]; % Do it for the other direction as well so that the two vertices pull each other.
end

% Generate an extraBonds matrix to make the tubules more rigid
for k = 1:numTubules
    for u = 1:tubuleHeight+1
        h = (k-1)*tubuleHeight + (k-1)*junctionSideLength + u; %Find the height of the ring of vertices we wan't to connect
        for i = 1:m
            q = [i + floor(m/2), h]; %Find the point that the vertex (i,h) is connecting to.
            q = makePeriodic(q, m);
            extraBonds{i,h} = q; %update it in the extraBonds matrix
            if mod(m,2) ~= 0 % If the circumference is odd
                r = q + [1,0]; % Then connect to a second point.
                r = makePeriodic(r, m);
                extraBonds{i,h} = [extraBonds{i,h}; r]; %update the matrix
            end
        end
    end
end

% Generate a matrix of additional bonds to push the junctions outwards
junctionExtraBonds = cell(1,1);
for k = 1:numTubules
    %Define the corners
    a1 = [c+(k-1)*x, k*tubuleHeight + (k-1)*junctionSideLength]; %define the reference corner a1
    a1 = a1 + [0,1]; %offset the height to match MATLAB's indexing
    a1 = makePeriodic(a1, m);
    a2 = a1 + [junctionSideLength, 0]; %Define the first 5-fold corner
    a2 = makePeriodic(a2, m); %Implement the periodic boundary
    a3 = a1 + [junctionSideLength, junctionSideLength]; 
    a3 = makePeriodic(a3, m);
    a4 = a1 + [0, junctionSideLength]; %Define the second 5-fold corner
    a4 = makePeriodic(a4, m);

    for j = 1:junctionSideLength-1
        for h = 1:junctionSideLength-1
            p = a1+[j,h]; %Find the point that we are pushing out.
            p = makePeriodic(p, m);
            junctionExtraBonds{p(1), p(2)} = [a1; a3]; %update it in the extraBonds matrix
        end
    end
end

end