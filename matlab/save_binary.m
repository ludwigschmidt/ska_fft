filename = '../../data/uv_8192_5_binary.dat';
to_store = complex_image;

file = fopen(filename, 'w');
fwrite(file, real(to_store), 'single');
fwrite(file, imag(to_store), 'single');
fclose(file);
