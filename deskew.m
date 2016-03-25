function [cbuf] = deskew(cbuf, nc, fdot, dr, c, pi,t_off)
%
% Filename:         deskew.m
%
% Removes a quadratic phase function from the deramp process
%
% cbuf = range data for ith pulse of phase history
% nc = range compression length
% fdot = chirp rate
% dr = range compressed sample spacing
% c = speed of light
% pi = value of pi
%==========================================================================
%==========================================================================
iswitch = 0;
r = zeros(nc,1);
p1 = zeros(nc,1);
p2 = zeros(nc,1);
p = zeros(nc,1);

s = (4*pi*fdot/c)/c;
drr = c*t_off/2;        % one way range difference between desired and actual
                        % mocomp point
if iswitch ==0
    s2 = 0;
else
    s2 = 8.0*pi*fdot*drr/c/c;   %ph_stab is used
end

i = 1:nc;
    r = dr*(i'-1-nc/2);      % dr = range sample spacing
    p1 = -s*r.^2;
    p2 = -s2*r;             % linear term needed only if using ph_stab
    p = p1+p2;
    cbuf(1:nc,1) = cbuf(1:nc,1).*(cos(p) + i*sin(p));
