% % % % % % % % % % % % % % % % % % % CALCULATE CENTROID OF THE INDIVIDUAL
% CELL
function [ x_com, y_com ] = CenterOfMass(I)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

global XMAX YMAX

CC = bwconncomp(I,8);   
 
z = regionprops(CC,'centroid','Area');  
 
 
 
 n = numel(z);
 
 

  M = 0;
 for i = 1:n
     M = M + z(i).Area;   
 end

 for i = 1:n
     theta_x(i) = (z(i).Centroid(1)/XMAX)*2*pi;
     theta_y(i) = (z(i).Centroid(2)/YMAX)*2*pi;
 end
%  z(2).Area
 
 for i = 1:n
     e1_x(i) = cos(theta_x(i));
     e2_x(i) = sin(theta_x(i));
     e1_y(i) = cos(theta_y(i));
     e2_y(i) = sin(theta_y(i));
 end
 
 
 
     mean_e1_x = 0;
     mean_e2_x = 0;
     mean_e1_y = 0;
     mean_e2_y = 0;
  
  for i = 1:n
     mean_e1_x = mean_e1_x +  (z(i).Area*e1_x(i))/M;
     mean_e2_x = mean_e2_x +  (z(i).Area*e2_x(i))/M;
     mean_e1_y = mean_e1_y +  (z(i).Area*e1_y(i))/M;
     mean_e2_y = mean_e2_y +  (z(i).Area*e2_y(i))/M;
  end
 
  mean_theta_x = atan2(-1*mean_e2_x,-1*mean_e1_x) + pi;
  mean_theta_y = atan2(-1*mean_e2_y,-1*mean_e1_y) + pi;
  
  x_com = (XMAX* mean_theta_x)/(2*pi); 
  y_com = (YMAX* mean_theta_y)/(2*pi); 


 
end

