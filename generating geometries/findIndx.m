function k = findIndx(i,h,m)
%FINDINDX Find the single integer index in the set of coordinates for (i,h)
%
% k = FINDINDX(i,h,m)
% 
% INTPUT:
%  i # First (horizontal) lattice parameter 
%  h # Second (vertical) lattice parameter 
%  m # Tubelet circumference, i.e., lattice periodicity
%
% OUTPUT:
%  k # Integer index for a vertex in the coordinates matrix

k = (h-1)*m + i;

end