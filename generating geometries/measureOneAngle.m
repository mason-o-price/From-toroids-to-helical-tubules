function [bindingAngle, N_A, N_B] = measureOneAngle(triangle_A, triangle_B, coordinates, m)
%MEASUREONEANGLE Measure the binding angle on a given edge
%
% [bindingAngle, N_A, N_B] = measureOneAngle(triangle_A, triangle_B, coordinates, m)
%
% INPUT:
%  triangle_A # Frist triangle containing the edge
%  triangle_B # Second triangle containing the edge
%  coordinates # Coordinates of the vertices
%  m # Tubelet circumference (T)
%
% OUTPUT:
%  bindingAngle # Angle (deg)
%  N_A % unit face normal on triangle A
%  N_B % unit face normal on triangle B


% Format: row = vertex, column 1 = h-coordinate of vertex, column 2 =
% k-coordinate of vertex

% Consider triangle A:
v1_A = triangle_A(1,:); % find vertex 1
v2_A = triangle_A(2,:); % vertex 2
v3_A = triangle_A(3,:); % vertex 3

% Now consider triangle B:
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
e3_A = v3_A_XYZ - v2_A_XYZ;

% Determine the (counter-clockwise) oriented edge vectors on triangle B
e1_B = v1_B_XYZ - v3_B_XYZ;
e2_B = v2_B_XYZ - v1_B_XYZ;
e3_B = v3_B_XYZ - v2_B_XYZ;

% Determine the unit face normal vectors:
N_A = cross(e1_A, e2_A)/norm(cross(e1_A, e2_A));
N_B = cross(e1_B, e2_B)/norm(cross(e1_B, e2_B));

center_A = (v1_A_XYZ + v2_A_XYZ + v3_A_XYZ)/3; % Find the center of triangle A
center_B = (v1_B_XYZ + v2_B_XYZ + v3_B_XYZ)/3; % Find the center of triangle B

d1 = norm(center_A - center_B);

offset_A = center_A + 0.2*N_A; % follow the normal vector away from the center of triangle A
offset_B = center_B + 0.2*N_B; % follow the normal vector away from the center of triangle B

d2 = norm(offset_A - offset_B);

if d2 < d1 % if the angle is negative (concave), then the offset vectors should become closer.
    sign = -1; % set the angle to be negative.
else %otherwise, assume the angle is positive.
    sign = 1;
end

% Calcualte the binding angle
bindingAngle = 180/pi*sign*abs(atan2(norm(cross(N_A, N_B)), dot(N_A, N_B)));
end