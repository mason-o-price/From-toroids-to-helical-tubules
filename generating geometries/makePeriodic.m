function newVertex = makePeriodic(vertex,m)
%MAKEPERIODIC Make the coordinates of a vertex periodic w.r.t. a given m
%
% newVertex = makePeriodic(vertex,m)
%
% INPUT:
%  vertex # original vertex coordinates [i,h]
%  m # lattice periodicity
%
% OUTPUT:
%  newVertex # updated vertex coordinates, with i between 1 and m.

%Implement the periodic boundary; taking 0 to be m.
newVertex = vertex;
newVertex(1) = mod(vertex(1)-1, m) + 1;
end