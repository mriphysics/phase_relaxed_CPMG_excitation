%% function to Gaussian convolve target magnetization in 3D
function P = GaussConv3D(P,w,m,X,Y,Z)

%%% don't always need m - can supply empty input
if numel(m)==0
    m=ones(size(P));
end

% Gaussian filter
R2 = X.^2+Y.^2+Z.^2;
Pg = exp(-R2/(2*w^2));

%%% depends on whether 3D or 4D (if using multi-freq.)
switch ndims(P)
    case {2,3}         %%% also works for 2D
        [mi,mii]=min(R2(:));
        [i1,i2,i3]=ind2sub(size(X),mii);
        Pg = circshift(Pg,[-i1 -i2 -i3]);
        Pgf = fftn(Pg);
        Pf = fftn(P);
        Pf = Pf .* Pgf;
        P = abs(ifftn(Pf));%<-- magnitude
        P=double(P/max(P(:)));
        % enforce mask
        P=P.*m;
    case 4

        R2=R2(:,:,:,1);
        [mi,mii]=min(R2(:));
        [i1,i2,i3]=ind2sub(size(X,1)*[1 1 1],mii);
        Pg = circshift(Pg,[-i1 -i2 -i3 0]);
        Pgf = ifft(fftn(Pg),[],4);
        Pf = ifft(fftn(P),[],4);
        Pf = Pf .* Pgf;
        P = abs(fft(ifftn(Pf),[],4));%<-- magnitude
        P=double(P/max(P(:)));
        % enforce mask
        P=P.*m;
end
end