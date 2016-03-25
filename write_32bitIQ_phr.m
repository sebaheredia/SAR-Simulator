function [] = write_32bitIQ_phr(im)
%
% Filename:         write_32bitIQ.m
%
% im = write_32bitIQ(im)
% Writes a raw signed 32-bit complex phase history file. 
% Data is written by column with alternating real and 
% imaginary 4-bit signed values
%
%INPUT:
% im        complex phase history to be written to file
%
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

name = input(['Input complex phase history filename to save: ']);
disp(['Choose path to save complex image file']);
path = uigetdir;
fname = [path '\' name]; % file to write to

type = 'float32'; 
img = im(:);
gg = length(img);
im_odd = real(img);
im_even = imag(img);

im_f(:,1) = zeros(2*gg,1);
im_f(1:2:end,1) = im_odd;
im_f(2:2:end,1) = im_even;

%Write the interlaced image to file 
fp = fopen(fname,'w','ieee-be');
if fp < 0
    fname
    error('Cant open float file')
end;
disp(sprintf('Writing complex data to %s ',fname));
count = fwrite(fp,im_f,'float32');
fclose(fp); 

