filename = '../../data/uv_16384_3_matlab_slabfft_result_binary.dat';
to_store = LL;

file = fopen(filename, 'w');
fwrite(file, real(to_store), 'single');
fwrite(file, imag(to_store), 'single');
fclose(file);
