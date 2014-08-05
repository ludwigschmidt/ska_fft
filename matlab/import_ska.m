%
% import_ska.m
%
% CJ Skipper, 2014-07-03
%
% import data from a CASA binary image file, and reformat it into a n-by-n
% array suitable for use in MATLAB. the CASA file format consists of an n-by-n
% array of complex numbers (4-byte float for real part and another 4-byte float
% for imaginary part), but is structured so that the data are stored in 32x32
% pixel blocks (from bottom-left to bottom-right, then working upwards).
%

% load datafile into a 1d array.
%fid = fopen( 'uv2048.dat' );
fid = fopen( 'uv16384_50.dat' );
raw_data = single(fread( fid, inf, 'float32' ));
fclose( fid );

% get the number of elements in the array.
num_elements = size( raw_data, 1 )

% calculate the number of pixels along each axis.
pixels_per_side = fix( sqrt( num_elements / 2 ) )

% CASA stores images in 32x32 pixel blocks, so calculate the number of blocks
% along each axis.
blocks_per_side = fix( pixels_per_side / 32 )

% reshape raw data into 2x? array (for real and imaginary components).
complex_data = reshape( raw_data, 2, [] );
clear raw_data;

fprintf( 'raw data reshaped.\n' )

% create a complex number array from the 2D array.
complex_array = complex( complex_data(1, :), complex_data(2, :) );
clear complex_data

fprintf( 'complex array loaded.\n' )

% organise the array into 4 dimensions: #3 and #4 to represent each 32x32 block
% and #1 and #2 to represent the position within the block. additionally,
% create a blank n-by-n image.
blocks = reshape( complex_array, 32, 32, blocks_per_side, blocks_per_side );
clear complex_array;
complex_image = zeros( pixels_per_side, pixels_per_side, 'single' );

fprintf( 'complex array reshaped.\n' )

% loop through each block, and then through the pixels in that block, and
% update the blank image with the correct values.
for k = 1:size( blocks, 3 )
	for l = 1:size( blocks, 4 )
		i = (k - 1) * 32;
		j = (l - 1) * 32;
		%for i = 1:32
			%for j = 1:32
				%complex_image( ((k - 1) * 32) + i, ((l - 1) * 32) + j ) = blocks( i, j, k, l );
		temp_matrix = squeeze( blocks(:, :, k, l) );
		complex_image([(i + 1):(i + 32)], [(j + 1):(j + 32)]) = temp_matrix([1:32],[1:32]);
			%end
		%end
	end
end
clear blocks;

fprintf( 'complex image built.\n' )

% finally, image should be transposed, and the y-axis flipped (since CASA
% stores images from the bottom up, with rows filled first, and MATLAB displays
% them from the top down, with columns filled first). the array 'complex_image'
% now contains the final data.
complex_image = flipud( permute( complex_image, [2,1] ) );

fprintf( 'complex image transposed and flipped.\n' )

% to display the image, I need to multiply each cell by some factor (since the
% values seem too small for the MATLAB colour scale) and then take the
% magnitude of each complex number.
figure; image( abs( complex_image * 10000000 ) );colormap('hot');

sky_image = fft2(complex_image);

figure; image(abs(sky_image * 40 ) );colormap('hot');

