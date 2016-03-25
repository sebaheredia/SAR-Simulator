function c = project(a,b)
%
% Filename:         bvec.m
%
% Projects vector a(*) into plane whose unit normal is b(*). Result is
% vector c(*).
%==========================================================================
%==========================================================================

s = dot(a,b);
for i=1:3
    c(i)=a(i)-s*b(i);
end 
	

