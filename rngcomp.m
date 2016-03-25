function [cbuf] = rngcomp(cbuf, nsamp, nc, inv, pi, n1, n2)
%
% Filename:         rngcomp.m
%
% Performs range compression by 1-D fft of length nc.  Data are
% symmetrically zero padded before fft.  Muxing is performed to shift DC
% term to center of array.
%
% cbuf = range data for ith pulse of phase history
% nsamp = number of samples in range dispersed pulse
% nc = range compression length
% inv = 0-- forward fft, 1--inverse fft
% pi = value of pi
%==========================================================================
%==========================================================================


%%%%% weight for sidelobe control %%%%%
% for i = 1:nsamp
%     wt = 0.5 + 0.5*cos(2*pi*(i-1-nsamp/2)/nsamp);   % Hann window
%     cbuf(i) = cbuf(i)*wt;
% end
wind=0;
ns = n2-n1+1;
if wind==1
    wt = 0.5 + 0.5*cos(2*pi*(i-n1-ns/2)/ns);   % Hann window
else
    wt = 1;
end

cbuf = cbuf(n1:n2,1).*wt;

%%%%% zero padding %%%%%
temp = zeros(nc,1)+i*zeros(nc,1);

temp(1:nsamp/2,1) = cbuf(nsamp/2+1:nsamp,1);
temp(nc-nsamp/2+1:nc,1) = cbuf(1:nsamp/2,1);   % swap data vector

cbuf(1:nc,1) = temp(1:nc,1);


%%%%% Muxing %%%%%
cbuf(2:2:nc) = -cbuf(2:2:nc);

%%%%% FFT %%%%%
%
if inv == 0
    cbuf = fft(cbuf,nc);
elseif inv ==1
    cbuf = ifft(cbuf,nc);
end 