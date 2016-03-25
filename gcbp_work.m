function image = cbp_work(n)
%
% Filename:         cbp_work.m
%
% This is the monostatic spotlight convolution backprojection image 
% formation algorithm. It allows for image formation in the ground 
% plane (in the local x-y-z system), the ground plane in the 
% xbar,ybar,zbar system or in the slant plane (in the local xhat-yhat-zhat system).
% Actually, the image can be formed in any plane defined by a ux,uy,uz
% coordinate system where the associated unit vectors are passed to
% the subroutines in place of ux,uy,uz.

% It assumes the input data is deramp-on-receive to the patch center
% for each transmit/receive pulse. The input complex data is
% filtered, range compressed via 1-d fft, deskew correction is performed,
% then each range-compressed and deskewed pulse is backprojected
% to a defined grid array. There are no limitations on where the grid
% is placed, the spacing between grid samples, or whether the grid is
% 2-D or 3-D. The user can write their own grid generator to use as
% necessary. This version assumes a 2-D grid layed out in the ground
% plane along the xbar-ybar axes, or in the slant plane along the 
% xhat-yhat slant plane axes.

% This program reads the transmitter path file created by
% the program pathgen.m and creates an array containing the unit vector
% triplets (x,y,z) pointing vectors from the imaged patch center to the
% monostatic SAR point for each transmit/receive pulse.
% Then, a specified transmit pulse is used to define the monostatic
% spheroid. A specified grid point (x,y,z) backprojects along the
% spheroid until it intersects somewhere along the monostatic line. The
% computed value of alpha specifies where along this line the 
% backprojection occurs in terms of +- meters from patch center. This
% value of alpha determines which samples in the corresponding range
% compressed and filtered pulse are used in the interpolation for
% subsequent backprojection.
%==========================================================================
%==========================================================================
disp('Specify phase history file to use in image formation. ');
phist = read_32bitIQ(400,408);

nsamp = n;
iswitch = 0;
c=3.0e8;		    % velocity of light in m/s

u_east = [1 0 0];   % initialize unit "east" "north" and "up" vectors
u_north = [0 1 0];
u_up = [0 0 1];

% Read transmitter path file, select starting and ending
% pulse numbers, and create new monostatic pointing unit vector array.
answt = 'n';%input(['Do you want Hann weighting of aperture data? (y,n) ']);
t_delay = 0;%input(['Specify additional time delay, t_delay = ']);	

disp('Specify transmitter path filename (for non-bistatic case use receive path file)');
[filename, pathname] = uigetfile('*');
tfile_path = [pathname filename];
tfile = dlmread(tfile_path);
[itcnt mm] = size(tfile);
xxt = tfile(:,2);
yyt = tfile(:,3);
zzt = tfile(:,4);
disp(['Number of triplets in transmitter file =',int2str(itcnt)]);

itstart = 1;%input(['Select starting pulse number in transmitter path file ']);
itend = 400;%input(['Select ending pulse number in transmitter path file ']);
ntp=itend-itstart+1;
disp(['Number of triplets selected in transmitter path file = ',int2str(ntp)]);
if(ntp > itcnt)
	error(['*** Too many pulses selected in transmitter file ***'])
end 

% Create new monostatic pointing vector array.
% The first triplet in the monostatic pointing vector array is to 
% correspond to the user-specified record in the phase history file.
clear pt;
pt(:,1) = xxt(itstart:itend,1);
pt(:,2) = yyt(itstart:itend,1);
pt(:,3) = zzt(itstart:itend,1);

% Compute various angles and vectors for diagnostic purposes.
% Specify ground plane unit normal in input coordinate system components
zbar = [0 0 1];
[xhat,yhat,zhat,xbar,ybar,theta_s,theta_g,slope,grazing,tilt] = angvec_dp(pt,ntp,zbar);
	
% Specify parameters that don't change with pulse number
disp('Note: The first paired tranmsitter/receiver position index that you selected must') 
disp('correspond to the proper phase history pulse record in the phase history file.')

nstart = 1;%input(['Specify the starting pulse number = ']);

% Specify a desired GRP offset vector.
disp('Note: the new GRP becomes the image grid origin.')
disp('The implied local North, East, and Up axes are unchanged.(Enter 0.,0.,0. for no offset.)')

dvgrp = [0 0 0];%input(['Specify the desired GRP offset vector [xg,yg,zg] = ']);
xg_off = dvgrp(1);yg_off = dvgrp(2);zg_off = dvgrp(3);

% Recompute pointing vectors and t_off based on a vector translation
% of the current GRP to a new GRP. Note: the new GRP becomes the origin
% of the image grid. The local North, East, and Up axes from which the
% original pointing vectors were computed remain unchanged with a
% shift of the GRP. This is consistent with a flat earth for translation
% purposes.
t_off = zeros(1,ntp);
rm=1e38;
for j=1:ntp
	rc2=0.0;            %mag squared of transmitter vector
	ptmag2=0.0;         %mag squared of new transmitter vector
    for i=1:3
		rc2 = rc2 + (pt(j,i)).^2;
		ptmag2 = ptmag2 + (pt(j,i)-dvgrp(i)).^2;;
    end
    
	rc=sqrt(rc2);
	ptmag=sqrt(ptmag2);
	if(ptmag < rm)		%get pulse and range of closest approach to the current GRP.
		rm=ptmag;
		jmin=j;
	end
    
	t_off(j)=t_off(j)-2.*(rc-ptmag)/c;   %overwrite t_off(*)
    for i=1:3
		pt(j,i)=pt(j,i)-dvgrp(i);	%overwrite vector pt(*)
	end 
end 
jmin=jmin-1+nstart;

disp('New pointing vectors computed.')
disp(['Pulse number of closest approach to GRP = ',int2str(jmin)]);
disp(['Range at closest approach to GRP = ',num2str(rm,6)]);
   
% Select ground-plane or slant-plane imagery options
disp('Note: For slant-plane imagery to be formed in the (xhat,yhat,zhat) coordinate system, enter "s".')
disp('For ground-plane imagery to be formed in the (east,north,up) coordinate system, enter "g"')
disp('For ground-plane imagery to be formed in the (xbar,ybar,zbar) coordinate system, enter "b"')
ansp = 'g';%input(['Enter image plane option: s, g, or b: '])

% Radar parameters
fo = 10e9;%input('Input radar center frequency (Hz): ') 
bwt = 4.015e8;%input('Input transmitter bandwidth (Hz): ') 
fdot = 1.0e12;%input('Input chirp rate (Hz/s): ') 
fs = 1.0e6;%input ('Input video a/d sampling rate after deramp (Hz): ')

%%%%% Data Parameters %%%%%
nsamp = 400;%input(['Specify number of input samples per pulse: ']);
disp(['There are ',int2str(nsamp),' samples per pulse.']);
n1 = 1;%input(['Specify starting sample number to process; n1 = ']);
n2 = 400;%input(['Specify ending sample number to process; n2 = ']);

%%%%% compute next power-of-two-larger size for range compression fft
nc = 2^(nextpow2(nsamp));
disp(['Selected FFT length for range compression = ',int2str(nc)])

%%%%% computed parameters %%%%%
t1=nsamp/2.0/fs;	%relative a/d start time (ahead of r/t time delay to patch center), t1=Tr/2. Effective receiver time, Tr=nsamp/fs.
bv = fs;                    % effective video bandwidth
dist = (c*bv/2)./(fdot);    % effective patch diamter (m)
disp(['Relative A/D start (ahead of patch delay time) = ',num2str(t1,6)])
dr = dist/nc;
disp(['The range compressed sample sapcing is ',num2str(dr),'.'])
disp(['Each range compressed pulse spans ',num2str(dist),'(meters).'])

%%%%% Output Data Parameters %%%%%
grid = [128 128];%input(['Specify output data grid array dimensions [nx ny nz]:  ']);
nx = grid(1); ny = grid(2); 

grid_s = [.1 .1];%input(['Specify output data grid sample spacings [dx dy dz]:  ']);
dx = grid_s(1); dy = grid_s(2); 

grid_off = [0 0];% %input(['Specify data grid center offset vector [xg_off yg_off zg_off]:  ']);
xg_off = grid_off(1); yg_off = grid_off(2); 

im = zeros(nx,ny)+i*zeros(nx,ny);
image = zeros(nx,ny)+i*zeros(nx,ny);
g_off = [xg_off yg_off];

% Write out unit vectors and aperture angles
disp(['Slant plane xhat = ',num2str([xhat(1),xhat(2),xhat(3)],4) ]);
disp(['Slant plane yhat = ',num2str([yhat(1),yhat(2),yhat(3)],4) ]);
disp(['Slant plane zhat = ',num2str([zhat(1),zhat(2),zhat(3)],4) ]);
disp(['Ground plane xbar = ',num2str([xbar(1),xbar(2),xbar(3)],4) ]);
disp(['Ground plane ybar = ',num2str([ybar(1),ybar(2),ybar(3)],4) ]);
disp(['Ground plane zbar = ',num2str([zbar(1),zbar(2),zbar(3)],4) ]);
disp(['Slant plane aperture angle (deg.) = ',num2str(theta_s,3) ]);
disp(['Ground plane aperture angle (deg.) = ',num2str(theta_g,3) ]);
disp(['Slope angle (deg.) = ',num2str(slope,3)]);
disp(['Grazing angle (deg.) = ',num2str(grazing,3) ]);
disp(['Tilt angle (deg.) = ',num2str(tilt,3) ]);

%%%%% define interval over which to find alpha(m) in pulse
rmin = -dist/2;
rmax = dist/2;

% Define one-sided sinc-function interpolator length. If lf=0, linear interpolation is performed.
lf=5;

% Begin processing pulses
cbuf = zeros(nc,1)+i*zeros(nc,1);
ict=0;
for npulse = nstart:ntp+nstart-1;%start:nstart+100
    ict=ict+1;
% Specify parameters that change for each pulse
    if(answt =='y')	    %Hann Weight	
		wt=.5+.5*cos(2.*pi*(npulse-nstart-ntp/2)/ntp);   %Hann weighting
		wt=2.*wt;	%keeps area under weight same as uniform
	else
		wt=1.0;	%uniform weight
	end 

% Form image in coordinate system defined by unit vectors ux,uy,uz. (ux,uy,uz)=(u_east,u_north,u_up), 
% (xhat,yhat,zhat), or (xbar,ybar,zbar)

    if(ansp == 's')     %convert pointing vectors
        ptnew = convert(pt(npulse-nstart+1,:),xhat,yhat,zhat);
	elseif(ansp == 'b')
		ptnew = convert(pt(npulse-nstart+1,:),xbar,ybar,zbar);
    elseif (ansp == 'g')	%default, if wrong answers to image plane selection
		ptnew = convert(pt(npulse-nstart+1,:),u_east,u_north,u_up);	
	end
    
    ptmag2 = zeros(1,3);    %mag squared of transmitter vector
    ptmag2 = sum(ptmag2 + ptnew.^2);
    ptmag = sqrt(ptmag2);
    rc = ptmag;	%one-way range to desired mocomp point
    ro = rc - c*t_off(npulse-nstart+1)/2.0;     %one-way range to actual mocomp point

% Read input data
    cbuf = phist(:,npulse);              % get pulse data from phase history
     
% Keep samples from n1 through n2
    cbuf(1:n1-1,:) = 0 + i*0;
    cbuf(n2+1:nsamp,1) = 0 + i*0;
   
    %%% Filter for subsequent backprojection
    [cbuf] = filter_s(cbuf, nsamp, fo, bwt, c);       % filter pulse data
    %%%%% range compress%%%%%
    inv = 1;                        % inverse fft
    [cbuf] = rngcomp(cbuf, nsamp, nc, inv, pi, n1, n2);
    %%%%% deskew %%%%%
    [cbuf] = deskew(cbuf, nc, fdot, dr, c, pi,t_off(npulse - nstart+1));

    im = gcbp_sub(im,cbuf,ptnew,g_off,dx,dy,wt,fo,ro,c,pi,rmin,rmax,dr,nc,lf);
    image = image+im;
 
    if mod(npulse,10)==0
        disp(['completed pulse',int2str(npulse)])
    else
    end
end 	%end pulse loop

[mm nn] = size(image);
disp(['The size of the computed complex imagery is = ',num2str([mm nn])])
write_32bitIQ(image);            % Save phase history

choice = input(['Do you want to disply imagery? (1--YES, 0--NO) ']);
if choice ==1
    imdisp(image',50,1);         % setup for correct axis when using imagesc
else 
end