function [ph_r,n] = sar_phgen
%
% Filename:         sar_phgen.m
%
% This function reads the transmitter/reciever path files created by
% pathgen.m and creates an ascii output file containing the (x,y,z) 
% pointing vectors from the imaged patch center to the bistatic sar
% point. This file is to be used as the pointing vector file in the
% subsequent processing and image formation.  This file also creates an 
% n x m complex file containing the phase history of multiple
% targets in spotlight radar format. The file writes pulses, each
% containing n complex fast time samples. There are m fixed length records 
% per file. The size of the file (n and m) depend on the various radar 
% parameters input to the program.
%
% The following is a list of the important variables used in the program:
% amp(i) = reflectance amplitude of the ith target
% dx(i) = change in x-position of the target per radar pulse (m)
% dy(i) = change in y-position of target per radar pulse (m)
% dz(i) = change in z-position of target per radar pulse (m)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% radar parameters (inputs) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fo = radar center frequency (Hz) 
% b = transmitted bandwidth (Hz) 
% fdot = chirp rate (Hz/s) 
% prf = pulse rep rate (Hz)
% fs = video a/d sampling rate after deramp (Hz)
% d = maximum slant plane patch diameter (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% computed parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmax = maximum valid range swath time interval after deramp
% t = transmit pulse length
% bx,by,bz = vector pointing from patch center back to bistatic radar point
% bsc = cos(beta/2) where beta = instantaneous bistatic angle
%==========================================================================
%==========================================================================
disp('')
disp('*****************************************************')
disp('*Bistatic Spotlight Radar Synthetic Target Generator*')
disp('******** with moving targets incorporated ***********')
disp('*****************************************************')
disp('')

CEE = 3.0e8; 
JAY = sqrt(-1);

disp('Select configuration file');
[filename, pathname] = uigetfile('*');
parfile = [pathname filename];
fid_par = fopen(parfile,'r');
fo = fscanf(fid_par,'%f',1); 
b = fscanf(fid_par,'%f',1);
fdot = fscanf(fid_par,'%f',1);
prf = fscanf(fid_par,'%f',1);
fs = fscanf(fid_par,'%f',1);
d = fscanf(fid_par,'%f',1);
signoise = fscanf(fid_par,'%f',1);
txfile = fscanf(fid_par,'%s',1);
rxfile = fscanf(fid_par,'%s',1);
tgtfile = fscanf(fid_par,'%s',1);
pointfile = fscanf(fid_par,'%s',1);
fclose(fid_par);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Input Parameters and Computed Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tp = 2*pi; c = 3.0e8;       % 2pi;velocity of light
beta = 0;                   %angel between transmit and recieve
t = b/fdot;                 % transmit pulse length
vb = 2.0*d*fdot/c;          % maximum video bandwidth after deramp(2*highest frequency)
tmax = t - 2.0*d/c;         % maximum valid pulse time interval

disp(['Maximum transmit pulse length is ' num2str(t,6)])
disp(['Maximum video bandwidth after deramp (min a/d rate) is ' num2str(vb,6)])
disp(['Maximum valid pulse time interval ' num2str(tmax,6)])
n = fix(fs*tmax);                % number of samples per pulse
disp(['Number of samples per pulse is ' int2str(n)])
fmul = 0;%input(['Specify fractional sample jitter multiplier: (values of 0-1) '])
iseed = 574904027;              % seed for random number generator
beff = fdot*tmax;           % effective bandwidth
q = fo/beff;                % Q factor (center freq/effective bandwidth)
disp(['Effective transmitted bandwidth is ' num2str(beff,6)])
disp(['System Q is ' num2str(q,6)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Get info from path files %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_tx = fopen(txfile,'r');
txin = fscanf(fid_tx,'%g %g %g %g %g',[5 inf]);
txin = txin';
xxt = txin(:,2); 
yyt = txin(:,3);
zzt = txin(:,4);
[itcnt mm] = size(txin);
fclose(fid_tx);

fid_rx = fopen(rxfile,'r');
rxin = fscanf(fid_rx,'%g %g %g %g %g',[5 inf]);
rxin = rxin';
xxr = rxin(:,2); 
yyr = rxin(:,3);
zzr = rxin(:,4);
[ircnt mm] = size(rxin);
fclose(fid_rx);
if itcnt ~= ircnt
        error('unequal length of transmitter and receiver path files')
end

bfile = 'bout.txt';%input('Specify filename for new bistatic pointing vector file ')
ntp = itcnt;
nrp = ircnt;
for i = 1:ntp
    xt = xxt(i);
    yt = yyt(i);
    zt = zzt(i);

    xr = xxr(i);
    yr = yyr(i);
    zr = zzr(i);

    pt = [xt yt zt];            
    pr = [xr yr zr];  
    [pb, bsc] = bvec(pt, pr, beta);     % computes bistatic pointing vectors
end
    file_name_ascii = [bfile '.txt'];
    save(file_name_ascii,'pb','-ascii');  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% read in target parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid_tgt = fopen(tgtfile,'r');
tgt_in = fscanf(fid_tgt,'%f %f %f %f %f %f %f %f',[8 inf]);
[jnk ntgt] = size(tgt_in);
%ntar = 1;
amp = tgt_in(1,:) + JAY*tgt_in(2,:);
xtgt = tgt_in(3,:);
ytgt = tgt_in(4,:);
ztgt = tgt_in(5,:);
dx = tgt_in(6,:);
dy = tgt_in(7,:);
dz = tgt_in(8,:);
fclose(fid_tgt);
[jnk ntar] = size(xtgt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Phase History Simulation%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% range computation %%%%%

ph_r = zeros(n,ntp) + i*zeros(n,ntp);
m = ntp;                % number of transmit pulses
for j = 1:m
    
    sbuf = zeros(1,n) + i*zeros(1,n);
        
    %%%%% calculate parameters that change each pulse %%%%%
    x = xxt(j);
    y = yyt(j);
    z = zzt(j);
    ptmag = sqrt(x.^2 + y.^2 + z.^2);       % inst. slant range to patch center from tx
    %
    x = xxr(j);
    y = yyr(j);
    z = zzr(j);
    prmag = sqrt(x.^2 + y.^2 + z.^2);       % inst. slant range to patch center from rx
    %
    ro = ptmag+prmag;           % total range to patch center
    ro2=ro*ro;
    to = -t/2 +(ro+d)/c;        % Video sample start time
    %%%%% Loop for targets %%%%%
    for i = 1:ntar
        xto = xtgt(i);          % initial positions of target
        yto = ytgt(i);
        zto = ztgt(i);
        xxx = xto + (j-1)*dx(i);    % positions due to target motion
        yyy = yto + (j-1)*dy(i);
        zzz = zto + (j-1)*dz(i);
        %
        rt2 = (xxt(j) - xxx).^2 + (yyt(j) - yyy).^2 + (zzt(j) -zzz).^2; %slant range to target
        rt = sqrt(abs(rt2));
        rr2 = (xxr(j) - xxx).^2 + (yyr(j) - yyy).^2 + (zzr(j) -zzz).^2;
        rr = sqrt(abs(rr2));
        r = rt+rr;                  % total range to target
        r2 = r*r;
        phi_1 = 2*pi*fo*(r-ro)/c;
        phi_2 = 2*pi*fdot*(r2-ro2)/(c*c);
        %
        %%%%% compute phase for target i along pulse j %%%%%
        k = 1:n;             % number of samples per pulse
        tt = to+(k-1)*tmax/n;       % time position along pulse j
        phi = mod(-2*pi*fdot*(r-ro)*tt/c - phi_1 + phi_2,tp);
        cp = cos(phi);
        sp = sin(phi);
        sbuf = sbuf+ amp(i).*complex(cp,sp); 

    end  
    ph_r(:,j) = sbuf'; 
end

[mm nn] = size(ph_r);
ph_r = ph_r ;%+ signoise*(randn(mm,nn)+j*randn(mm,nn));       % add WGN

disp(['The size of the phase history data is = ',num2str([mm nn])])
write_32bitIQ_phr(ph_r);            % Save phase history

%% compute image using 2D FFT and display
choice = input(['Do you want to compute imagery using 2-D FFT and dispay? (1--YES, 0--NO) ']);

if choice ==1
    load c_map.dat;
    imagesc(real(ph_r));colormap(c_map);
    sar_image = fftshift(abs(ifft2(ph_r,256,256)).^2);
    dB_range = input('Select top dynamic range in dB to plot: ');
    log_sar = 10*log10(sar_image);
    log_sar = log_sar - max(max(log_sar));
    vec = find(log_sar < (-dB_range));
    log_sar(vec) = -dB_range*ones(size(vec));
    figure;imagesc(log_sar);colormap(c_map);
else 
end