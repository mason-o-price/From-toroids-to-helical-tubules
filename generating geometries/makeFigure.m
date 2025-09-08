function makeFigure(i)
%MAKEFIGURE Prepare a figure for 3D plotting
%
% MAKEFIGURE(i)
%
% INPUT:
%  i # figure index

clf(figure(i))
figure(i);
hold on; %make sure the plots don't overwrite each other
axis off; %turn off the grid
set(gcf,'Color','white'); %make the background white
set(gcf, 'renderer','painters'); %make sure we can output vector file figures
axis equal %set the axes to be equal size
rotate3d(figure(i),'on');

end