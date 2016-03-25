function [cbuf] = filter_s(cbuf, nsamp, fo, bwt, c)
%
% Filename:         filter.m
%
% Function to apply the frequency-domain "filtering" on the range-dispersed
% pulse for subsequent backprojection.

% cbuf = range data for ith pulse of phase history
% nsamp = number of samples in range dispersed pulse
% fo = center carrier frequency
% bwt = transmitter bandwidth
% c = speed of light
%==========================================================================
%==========================================================================
x = zeros(nsamp,1);
rho = zeros(nsamp,1);

i = 1:nsamp;
x = (i'-1-nsamp/2)/(nsamp/2);           % number between -1 and 1
rho = x*bwt/2;                          % fraction of bandwidth
cbuf(1:nsamp,1) = cbuf(1:nsamp,1).*abs(rho+fo)*2/c;      % selects out range samples within transmitter bandwidth
cbuf(1:nsamp,1) = cbuf(1:nsamp,1).*abs(rho+fo)/fo;      % normalized


