% Script - SAR Simulator
% Filename:         sar_simulator.m
% Author:           drohm
% This script is used to initialize and call the functions for simulating a
% synthetic aperture radar collection.  
%==========================================================================
%==========================================================================
clear all;close all;clc

%% Input values for radar platform path generator
npass = 1;  %input(['Input the number of passes for this volumetric collection: ']);
elevsp = 1; %input(['Input elevation spacing between passes: '])
    
    vx = 500;           % x-velocity of platform (m/sec)
    vy = 0;             % y-velocity of platform (m/sec)
    vz = 0;             % z-velocity of platform (m/sec)
    prf = 163.5;        % pulse (xmit or recieve) repitition frequency (Hz).
    tmax = 2.5;         % maximum time extent of path (sec)
    xo = -613;          % initial x-position of path (m) 
    yo = -25980;        % initial y-position of path (m)
    zo = 15000;         % initial z-position of path (m)

for i = 1:npass
    zo = zo + (i-1)*elevsp;
    [pathg] = pathgen(vx,vy,vz,prf,tmax,xo,yo,zo);      % Funtion call 

    %% Phase History Generator
    [phr,n] = phgen;           % Function call - phase history generator
end
%% Convolution Back Projection
for i = 1:npass
    image = gcbp_work(n);        % convolution backprojection image formation
end
