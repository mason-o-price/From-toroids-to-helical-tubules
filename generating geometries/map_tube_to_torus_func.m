function new_coordinates = map_tube_to_torus_func(T,L,D,R,numSegments,tube_radius,coordinates,rate_offset)
%MAP_TUBE_TO_TORUS_FUNC Map vertices on a cylinder to a torus
% 
% new_coordinates = map_tube_to_torus_func(c,L,d,x,numSegments,tube_radius,coordinates,rate_offset)
% 
% INPUT:
%  T # Tubelet circumference
%  L # Tubelet length
%  D # Junction side-length
%  R # Rolling parameter
%  numSegments # Number of tubelet segments N
%  tube_radius # Radius of the tubelet
%  coordinates # Matrix of coordinates
%  rate_offset # Offset twist to tune the location of vertices on the torus
%
% OUTPUT:
%  new_coordinates # Matrix of coordinates on the torus

%% Preparation
% Rename variables
c = T;
d = D;
x = R;

% Find the total height, (toroid's circumference)
totalHeight = numSegments*(L + d) + 1;
% total height of the initial tubule
H = (totalHeight-1)*sqrt(3)/2;

% torus parameters
resolution = 50;
[alpha,phi] = meshgrid(linspace(0,2*pi,resolution));
r2 = tube_radius;
R = H/(2*pi);
x_torus = (R + r2*cos(alpha)).*cos(phi);
y_torus = (R + r2*cos(alpha)).*sin(phi);
z_torus = r2*sin(alpha);
K = cos(alpha)./(r2*(R+r2*cos(alpha)));

% start with xyz coordinates and convert to cylindrical coordinates: h,k,theta,r,z
cylinder_coordinates = coordinates;
cylinder_coordinates(:,3) = atan2(coordinates(:,4), coordinates(:,3)); % atan2 returns 
cylinder_coordinates(:,4) = sqrt(coordinates(:,3).^2 + coordinates(:,4).^2);
cylinder_coordinates(:,5) = -coordinates(:,5);

% define parameters to adjust the map
side3_bevelAngles = [-38.9424, -27.6723, -21.6246, -17.7989, -15.1445, -13.189, -11.6859, -10.4933, -9.5233, -8.7185, -8.0398, -7.4596, -6.9579, -6.5197, -6.1335, -5.7907, -5.4842, -5.2087]; % starting with (3,0) and assuming achiral
axial_height = sqrt(3)/2; %sqrt(3)/2*cosd(side3_bevelAngles(c-2));
alpha_offset = L*2*pi/c +pi/4; %+pi/4;
%rate_offset = 0.5; %offset the rate of twist w.r.t z

% map the cylindrical coordinates to the torus, (i.e. coordinates u,v)
torus_coordinates = cylinder_coordinates(:,1:2);
torus_coordinates(:,3) = cylinder_coordinates(:,5)/R;
torus_coordinates(:,4) = cylinder_coordinates(:,3) + cylinder_coordinates(:,5)/(axial_height*(L+d))*(1+rate_offset)*pi*x/c+alpha_offset;

% now use the parametrization of a torus to find the xyz coordinates of the vertices
new_coordinates = torus_coordinates;
new_coordinates(:,3) = (R+r2*cos(torus_coordinates(:,4))).*cos(torus_coordinates(:,3));
new_coordinates(:,4) = (R+r2*cos(torus_coordinates(:,4))).*sin(torus_coordinates(:,3));
new_coordinates(:,5) = r2*sin(torus_coordinates(:,4));

end