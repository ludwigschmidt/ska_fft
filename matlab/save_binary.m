filename = 'uv_2048_binary.dat';

file = fopen(filename, 'w');
fwrite(file, real(complex_image), 'single');
fwrite(file, imag(complex_image), 'single');
fclose(file);