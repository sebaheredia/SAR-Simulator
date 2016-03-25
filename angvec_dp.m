function [xhat,yhat,zhat,xbar,ybar,theta_s,theta_g,slope,grazing,tilt] = angvec(pt,ntp,zbar)
%
% Filename:         angvec.m
% 
% Routine to compute slant plane and ground plane unit vectors from
% selected pointing vectors. The slant and ground plane aperture
% angles, slope, grazing, and tilt angles are also computed.
% zbar(*) is input as the ground plane unit normal in input
% coordinate system components. The slant plane will be defined by the
% pointing vectors at 20% and 80% through the selected aperture range.
% The slant plane mid-aperture pointing vector will be used to define
% yhat.
%==========================================================================
%==========================================================================

%xhat = zeros(1,3);yhat = zeros(1,3);zhat = zeros(1,3);	%slant plane unit vectors
%xbar = zeros(1,3);ybar = zeros(1,3);zbar = zeros(1,3);	%ground plane unit vectors
%temp1 = zeros(1,3);temp2 = zeros(1,3);temp3 = zeros(1,3);temp4 = zeros(1,3);
mid=fix(ntp/2+ .5);	%mid-aperture index
m20=fix(ntp*0.2 + .5);	%20% aperture index
m80=fix(ntp*0.8+ .5);	%80% aperture index
%
%pt = pt';
for i=1:3
    temp3(i)=pt(mid,i);
	temp1(i)=pt(m20,i);
	temp2(i)=pt(m80,i);
end 
temp3 = temp3./norm(temp3);
temp1 = temp1./norm(temp1);
temp2 = temp2./norm(temp2);
zhat = cross(temp1,temp2);  %slant plane unit normal
zhat = zhat./norm(zhat);
% Project mid-aperture unit vector into slant plane to define yhat
yhat = project(temp3,zhat);
yhat = yhat./norm(yhat);
%
% Make sure slant plane normal has a component in the same direction as zbar.
x = dot(zhat,zbar);
if(x<0e0)
	for i=1:3
        zhat(i) = -zhat(i);
    end
else
end
%%Get ground plane ybar by projecting yhat into ground plane
ybar = project(yhat,zbar);
ybar = ybar./norm(ybar);	    %ground plane ybar
%
%Compute xhat and xbar
xbar = cross(ybar,zbar);
xbar = xbar./norm(xbar);  	%ground plane xbar
xhat = cross(yhat,zhat);
xhat = xhat./norm(xhat);	    %slant plane xhat
%
%	Compute slant and ground plane aperture angles
for i=1:3
	temp1(i)=pt(1,i);	%first pointing vector
	temp2(i)=pt(ntp,i);	%last pointing vector
end 
temp1 = temp1./norm(temp1);
temp2 = temp2./norm(temp2);
%
% Project first and last pointing vectors into slant plane to
% compute slant plane aperture angle.
temp3 = project(temp1,zhat);
temp3 = temp3./norm(temp3);	    %first slant plane vector
temp4 = project(temp2,zhat);
temp4 = temp4./norm(temp4);	    %last slant plane vector
%
cts = dot(temp3,temp4);
theta_s=(acos(cts))*180/pi;	        %slant plane  aperture angle (deg.)
%
%% Project first and last pointing vectors into ground plane to
%% compute ground plane aperture angle.
temp3 = project(temp1,zbar);
temp3 = temp3./norm(temp3);	    %first ground plane vector
temp4 = project(temp2,zbar);
temp4 = temp4./norm(temp4);        %last ground plane vector

ctg = dot(temp3,temp4);
theta_g = (acos(ctg))*180/pi;	    %ground plane  aperture angle (deg.)
%
%%Compute slope, grazing, and tilt angles.
cs = dot(zhat,zbar);
cg = dot(yhat,ybar);
ct = dot(xhat,xbar);
slope = (acos(cs))*180/pi;
grazing = (acos(cg))*180/pi;
tilt = (acos(ct))*180/pi;

for i = 1:3
    xhat(i) = (xhat(i));
    yhat(i) = (yhat(i));
    xbar(i) = (xbar(i));
    ybar(i) = (ybar(i));
end