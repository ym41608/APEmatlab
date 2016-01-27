function [trans_mat, ex_mat] = getTransAndExMatrix(pose,in_mat)

tx = pose(1);
ty = pose(2);
tz = pose(3);
rx = pose(4);
rz0 = pose(5);
rz1 = pose(6);

% Calculate rotation matrix R and translation vector t
%R = [cos(rz0) -sin(rz0) 0; sin(rz0) cos(rz0) 0; 0 0 1] * ...
%	 [1 0 0; 0 cos(rx) -sin(rx); 0 sin(rx) cos(rx)] * ...
%	 [cos(rz1) -sin(rz1) 0; sin(rz1) cos(rz1) 0; 0 0 1];
R = rotz(rad2deg(rz0))*rotx(rad2deg(rx))*rotz(rad2deg(rz1));
R(:,3) = -R(:,3); 		%%  y cross x
t = [tx ty tz]';

% Calculate extrinsic matrix
ex_mat = [R t];

% Calculate transformation matrix 
trans_mat = in_mat * [ex_mat; 0 0 0 1];



