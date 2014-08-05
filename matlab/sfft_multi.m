n=16384;
w=64;

projh = fft2(complex_image(floor((n-w)/2)+1:floor((n+w)/2),:));
projv = fft2(complex_image(:,floor((n-w)/2)+1:floor((n+w)/2)));

th = 0.3;

b = n/w;

Lh = abs(projh) > th;
Lv = abs(projv) > th;

[Xh Yh] = find(Lh > 0);
[Xv Yv] = find(Lv > 0);

pr = -3*b/2:1:-b/2;

peaks = ones(max(length(Xv),length(Xh))*b,2);
ind = 1;
for ii =1:1:length(Yh)
% NEED TO TAKE CARE OF BOUNDARIES


YvP = round((Yh(ii))/b)+1;
XvP = Xh(ii)*b+pr;

Xc =intersect(Xv(Yv==YvP),XvP);

lc = length(Xc);

candidates = [Xc.',Yh(ii)*ones(lc,1)];

peaks(ind:ind+lc-1,:) = candidates;
ind = ind + lc;
end

peaks = peaks(1:ind,:);

LL = zeros(n,n);
for ii=1:1:ind
    v = peaks(ii,1);
    h = peaks(ii,2);
    %LL(v,h) =0.5*abs(projh(round(v/b)+1,h))+0.5*abs(projv(v,round(h/b)+1));
    LL(v,h) =sqrt(abs(projh(round(v/b)+1,h)*projv(v,round(h/b)+1)));
end

figure(1);imagesc(abs(LL));colormap('hot');

