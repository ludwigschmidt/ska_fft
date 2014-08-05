n = 16384;
w = 64;
w2 = 1024;
%n = 2048;
%w = 8;
%w2 = 128;


% Middle Square
projs = fft2(complex_image(floor((n-w2)/2)+1:floor((n+w2)/2),floor((n-w2)/2)+1:floor((n+w2)/2)));

% Horizontal and Vertical Slices
projh  = fft2(complex_image(floor((n-w)/2)+1:floor((n+w)/2),:));
projv  = fft2(complex_image(:,floor((n-w)/2)+1:floor((n+w)/2)));


% Diagonal Slices
projt1 = zeros(w,n);
projt2 = zeros(w,n);
for ii=1:1:n
    projt1(:,ii) = complex_image(mod(ii-1-w/2:ii+w/2-2,n)+1,ii);
    projt2(:,ii) = complex_image(mod(-(ii-1-w/2:ii+w/2-2),n)+1,ii);
end
projd1 = fft2(projt1);
projd2 = fft2(projt2);

% Threshold (Set by trail and error)
th = 0.35;

% Number of frequencies per bin
b = n/w;
b2 = n/w2;

% Large Bins  in Mid square slice
Ls = abs(projs) > th;
[Xs Ys] = find(Ls > 0);

% Frequency range in the bin
pr = -3*b2/2:1:-b2/2;

% Number of large bins
P = length(Xs);

LL = zeros(n,n);
for ii=1:1:P
    
    x = Xs(ii);
    y = Ys(ii);

    % List frequencies in large bin
    xr = mod(x*b2 + pr-1,n) +1;
    yr = mod(y*b2 + pr-1,n) +1;
    
    % For each frequency in large bin, estimate it's value from the other binnings.
    for vr = 1:1:length(xr)
    for hr = 1:1:length(yr)
        v = xr(vr);
        h = yr(hr);

        As  = abs(projs(mod(round(v/b2),w2)+1,mod(round(h/b2),w2)+1));
        Ah  = abs(projh(mod(round(v/b),w)+1,h));
        Av  = abs(projv(v,mod(round(h/b),w)+1));
        Ad1 = abs(projd1(mod(round(v/b),w)+1,mod(v-1+h-1,n)+1));
        Ad2 = abs(projd2(mod(round(h/b),w)+1,mod(-v+1+h-1,n)+1));
       
        LL(v,h) = mean([As,Ah,Av,Ad1,Ad2]);
        
    end
         
    end
    
end

% Plot result
figure;imagesc(abs(LL));colormap('hot');

