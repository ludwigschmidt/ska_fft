function res = load_slabfft_result(filename, n)
    f = fopen(filename, 'r');
    [re, count] = fread(f, [n n], 'single');
    if count ~= n * n
        error('Read only %d real values.', count);
    end
    [im, count] = fread(f, [n n], 'single');
    if count ~= n * n
        error('Read only %d imaginary values.', count);
    end
    fclose(f);
    res = complex(re, im);
end

