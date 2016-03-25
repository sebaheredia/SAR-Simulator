function imdisp(x,dB,mag)
%
% Filename:         imdisp.m
%
% function call to display image x on logorithmic scale.  Also
% computes and dispays spectum.
% mag = 1 if x is complex valued
% mag = 0 if x is absolute valued
%
%==========================================================================
%==========================================================================
load c_map.dat;
dB_range = dB;
if mag ==1
    im_dB = 10*log10(abs(x).^2);
elseif mag ==0
    im_dB = 10*log10(x);
end
im_dB = im_dB - max(max(im_dB));
vec = find(im_dB < (-dB_range));
im_dB(vec) = -dB_range*ones(size(vec));
figure;imagesc(im_dB);colormap(c_map);
axis tight
daspect([1 1 1])
