function [pb, bsc] = bvec(pt, pr, beta)
%
% Filename:         bvec.m
%
% Function bvec.m computes the bistatci pointing vector, pb, 
% which bisects the angle between the transmitter pointing 
% vector, pt, and the reciever pointing vector, pr, and 
% lies in the plane formed by pt and pr.  The magnitude of pb
% is the average of the magnitude of pt and pr.
%
% pt = transmit pointing vector
% pr = reciever pointing vector
% beta = angle between pt and pr
%==========================================================================
%==========================================================================
eps = 1.0e-6;
for i = 1:3
    temp1(i) = pt(i);
    temp2(i) = pr(i);
end
prmag = abs(pr).^2;
ptmag = abs(pt).^2;
temp1 = norm(pt);
temp2 = norm(pr);
pb = temp1+temp2;
pbmag = abs(pb).^2;

if pbmag<=eps
        pb(i)==0.0;
        bsc==0.0
else 
end
pb = norm(pb);
bsc = dot(temp1,pb);
sc = (abs(pt) + abs(pr))/2.0;
pb = pb*sc;
%pb = (pt/abs(pt) +pr/abs(pr))*sc;
%bsc = cos(beta/2);