function im = read_32bitIQ(num_col,num_row)
%
% Filename:         read_32bitIQ.m
%
% mod 2/21/05 -  Modified to read arbitrary sized image chip from file
%
% im = read_32bitIQ(num_col, num_row)
% Reads a raw signed 32-bit complex image file. 
% Image has width num_col pixels and hieght num_row.  
% Assumes the data is written by column with alternating real and 
% imaginary 4-bit signed values.
% Each pixel is the sum of a real and complex value
%
% INPUT:
% num_col:		Number of columns in image 0 -> M
% num_row:	    Number of rows in image 0 -> N
%
% OUTPUT:
% im:		The complex 32-bit IQ image
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
[filename, pathname] = uigetfile('*','Select image file for display');
filename = [pathname filename];

if( num_col<0 ) 
  error('x indices out of range')
end
if( num_row<0 )
  error('y indices out of range')
end
Ncol = num_col;
disp(sprintf('  Reading from: %s',filename));
fp = fopen(filename,'r','ieee-be');
if fp == -1 
    error('Cant open file'); 
end
sp = 2*Ncol*Ncol;  
im = zeros(num_row,num_col);
%fseek(fp,'BOF',0);
fseek(fp,0,'bof');

data = fread(fp,sp,'float32');     
imc = data(1:2:end,1) +j*data(2:2:end,1);
im = reshape(imc,Ncol,Ncol);
fclose(fp); 

