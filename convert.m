function [ptnew] = convert(pt,ux,uy,uz)
%
% Filename:         convert.m
%   
% Routine to convert pointing vector defined in a local z,y,z
% basis to pointing vector defined in a local ux,uy,uz basis.
%
%==========================================================================
%==========================================================================
ptnew = zeros(1,3);
ptnew(1) = dot(pt,ux);  %Put ux component in ptnew(1)	
ptnew(2) = dot(pt,uy);  %Put uy component in ptnew(2)	
ptnew(3) = dot(pt,uz);  %Put uz component in ptnew(3)
	
%ptnew = ptnew./norm(ptnew);