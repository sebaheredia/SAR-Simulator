function [pathg] = pathgen(vx,vy,vz,prf,tmax,xo,yo,zo)
% Filename:         pathgen.m

% This function creates an ascii file containing triplets (x,y,z)
% describing a path in 3-dimensions. The number of triplets
% generated depends on the input parameters. This function is used to
% create a transmitter path and/or reciever path for use by the
% bistatic/monostatic synthetic target generator (sar_phgen.m). The
% following input variables are used in the function:

% vx = x-velocity of platform (m/sec)
% vy = y-velocity of platform (m/sec)
% vz = z-velocity of platform (m/sec)
% prf = pulse (xmit or recieve) repitition frequency (Hz).
% tmax = maximum time extent of path (sec)
% xo = initial x-position of path (m)
% yo = initial y-position of path (m)
% zo = initial z-position of path (m)

%==========================================================================
%==========================================================================
%
nmax = prf*tmax;    %maximum number of pulses transmitted/recieved
xmax = vx*tmax;     %distance in x,y,z dimensions platform transmitted/recieved
ymax = vy*tmax;         
zmax = vz*tmax;
dt = 1.0/prf;       % Period of pulses (time between pulses)
%
disp(['The total number of triplets (x,y,z) to be computed is ', int2str(nmax)])
disp(['The total distances traveled in x,y,z are [', int2str([xmax,ymax,zmax]),']'])
%
tjit = 0;
for i = 1:nmax
    count(i) = i;
    x(i) = xo+vx*dt*(i-1);
    y(i) = yo+vy*dt*(i-1);
    z(i) = zo+vz*dt*(i-1); 
    bnk(i) = 0.0;
    
end  
pathg = [count' x',y',z' bnk'];
%
%%%%% Draw Geometry %%%%%
sx = [0 xo./1000];sy = [0 yo./1000];sz = [0 zo./1000];
ex = [0 x(end)./1000];ey = [0 y(end)./1000];ez = [0 z(end)./1000];
plot3(sx',sy',sz','k');hold on;
plot3(ex,ey,ez,'k')
plot3(x./1000,y./1000,z./1000,'bo');box on;grid on
xlabel('x-dir (km)');ylabel('y-dir (km)');zlabel('Elevation (km)')
plot3(0,0,0,'ko');
hold off

%%%%% Save path file %%%%%
choice = input(['Do you want to save this flight path? (1--YES,0--NO) ']);
if choice ==1
    disp(['Choose a directory to save path file']);
    pathname = uigetdir;
    file_name = input(['Input file name for saved path file (ascii - use single quotes):  ']);
    file_name_ascii = [pathname '/' file_name '.txt']
    save(file_name_ascii,'pathg','-ascii');  
else
end

    
  