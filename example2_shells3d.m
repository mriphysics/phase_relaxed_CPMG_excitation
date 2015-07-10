% Example 3D shells trajectory phase relaxed 90° pulse.
% Shaihan Malik 10th July 2015. 

%% constant definitions
gamma_uT = 267.5221;        % units rad s^-1 uT^-1
gamma_mT = gamma_uT * 1e3;
dt = 6.4e-6;                % RF dwell time, 6.4us

flip = 90;                  % total flip angle, degrees
flip=flip*pi/180;

%%% Code depends on VERSE guided time optimal design (doi: 10.1002/mrm.21950)
%%% This is implemented and available as another git repo from @mriphysics. 
%%% See the readme for more details. 
addpath(genpath('../reverse-GIRF/'))

%%% Also library functions
addpath('lib')

%% Load in B0 and B1 maps, plus mask

% B0 is in units of uT. B1 field ('tx') is relative sensitivity
load b0_b1_maps_3d
[n1 n2 n3 Nc]=size(tx);           

% random downsample - in order to test on a computer without sufficient
% memory, can drop random voxels (this was NOT used in the publication;
% comment these two lines out to remove).
rr = rand(size(m));
m = m & (rr<0.3);

% Make list of coordinates, and index of those within mask
xlist = double(cat(2,X(:),Y(:),Z(:)));
idx = find(m);
Nidx=length(idx);



%% target excitation

%%% offset and size
xoff=-4.7e-3;yoff=-1e-3;zoff=-0.5e-3;
r0=0.03;

%%% Sphere
%P = ((X-xoff).^2 + (Y-yoff).^2 + (Z-zoff).^2)<r0^2;

%%% Cube
P = double(abs(X-xoff)<r0 & abs(Y-yoff)<r0 & abs(Z-zoff)<r0);

%%% apodize
P = GaussConv3D(P,6.5e-3,m,X,Y,Z);

%%% Now add pi/2 phase (CPMG) and scale to flip
P=P*flip*1i;


%% Load in GIRF and prepare distortion model

% GIRF
girf = load('GIRF_3T_London_20140729');

%%% FWD model
gcor =  @(x)(gradient_distort_GIRF(x,girf.ff,girf.Hw,dt,10));


%% Define k-space trajectory

%%% Generate each shell separately using Wong & Roos method
k = {};
K=[];
R = [30:50:300];
R=fliplr(R);
nr=length(R);

%%% each shell oriented differently
ra = pi/6*((1:nr)-1);

%%% surface area prportional to r^2
num = ceil(150*(R(:)/R(1)).^2);
min_num = 30;
num(num<min_num)=min_num;%<--- minimum value
disp(num')
for ii=1:nr
    k{ii} = R(ii)*gen_single_shell(num(ii),'rot',ra(ii)+pi*mod(ii,2),'slow',0.00005);
    K = cat(1,K,k{ii});
end

%%% explicitly set last value to zero
K = cat(1,K,zeros([1 3]));

%%% First gradient trajectory - this only uses gradient constraints (not RF)
%%% Define target trajectory using Wong & Roos method then make time
%%% optimal waveforms using Lustig method (10.1109/TMI.2008.922699)

Kcm = K / (2*pi*100);% k-space in 1/cm
Gmax = 40 / 10; % G/cm
Smax = 180 / 10; % G/cm/ms
[C,time,g,s,k, phi, sta, stb] = minTimeGradient(Kcm,0,0,0,Gmax,Smax,dt*1e3,[],[]);
fprintf(1,'Duration = %1.3f ms\n',time);

%%% check actual TRANSMIT kspace
G = 10*g; %<- convert back to mT/m
M=length(G);

%%% simulate the effect of Eddy currents 
[Gedd,kedd] = gcor(G);
    
t = (1:M)*dt;
klist=kedd;



%% reVERSE
%%% For convenience we also use MCGLS for inverting the linear problem:
%%% (http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=3&SC_ID=10&MP_ID=126)
%%% published here: doi:10.1016/j.apnum.2003.11.010
%%% Replace with any regularized linear inverse if mcgls not available.

%%% Pick fixed Wreal and lambda
Wreal_use = 0.25;
lambda_use = 100;

%%% Options for reVERSE
verseopt = reVERSE_init;
verseopt.b1_limit = 11 * 1e-2; % 11uT absolute limit (uT -> Gauss)
verseopt.bmax =     8 * 1e-2;  % 8uT target on each VERSE iteration (uT -> Gauss)
verseopt.os = 15; %<--- oversample factor

% initialize gradients and k-space
Kv = klist;
kv=k;
tv=t;
gv = g;
Mv=length(Kv);

% variables
bb = {};
PP = {};
bv = {};
Gv = {};

% 1st gradient is already computed
Gv{1}=gv*10;%<-- convert to mT/m
Nstop=5;

for ii = 1:Nstop
    %%% Regularized STA design
    fprintf(1,'iteration %d: define STA system matrix\n',ii);
    A = STA_system_matrix(xlist,Kv,tv,tx,b0,m,'loopcalc');
   
    fprintf(1,'iteration %d: split RE and IM\n',ii);
    %%% Now split real and imaginary parts
    Aw = [imag(A) real(A);real(A) -imag(A)];
    Pw = [imag(P(idx));real(P(idx))];
        
    %%% can add a weighting in here for the real/imag channels
    WW = (Wreal_use+1i)*double(m);
    ww = [imag(WW(idx));real(WW(idx))];
    Aw = Aw .* repmat(ww,[1 size(Aw,2)]);
    Pw = ww.*Pw;
    
    % Use best lambda from above.
    fprintf(1,'iteration %d: Solve LLS problem \n',ii);
    
    opt=struct('max_it',50);
    [b,brelres] = mcgls(Aw,Pw,lambda_use,opt);
    fprintf(1,'... done. Relative residual %1.3f, max B1 %1.2f uT\n',brelres(end)/norm(Pw(:)),1000*max(abs(b(:))));
    %%% now split
    b = complex(b(1:Mv*Nc),b(Mv*Nc+1:Mv*Nc*2));
    b = reshape(b,[Mv Nc]);
    clear Aw Pw

    % put into real space
    Psol = zeros(size(P));
    Psol(idx) = A * b(:);
    
    % save these
    bb{ii}=b;
    PP{ii}=Psol;

    % VERSE: Only do this if there's another iteration
    if (ii==Nstop)||(max(abs(b(:))*10)<verseopt.b1_limit)
        break
    end
    fprintf(1,'iteration %d: VERSE\n',ii);
    
    %%% Update optional inputs with RF pulse and Grads (both in CGS units)
    verseopt.b = b*10;   %<-- mT->Gauss
    verseopt.Gmod = sum(gv.^2,2).^0.5;

    [Cv,time,gv,sv,kv, phi, sta, stb,p_of_t,oo] = minTimeGradient_VERSE(kv,0,0,Gmax,Smax,dt*1e3,[],[],verseopt);

    %%% get bv from function
    bv{ii} = oo.bt/10;% back to mT

    %%% Re-make Tx k-space (this is for next STA run)
    Gv{ii+1} = 10*gv; %<- convert back to mT/m
    Mv=length(Gv{ii+1});
    [~,kedd] = gcor(Gv{ii+1});
    tv = (1:Mv)*dt;
    Kv = kedd;
    
    figure(1);
    nr=2;nc=2;sl=16;
    subplot(nr,nc,1)
    imagesc(abs([imag(Psol(:,:,sl)) real(Psol(:,:,sl))]))
    
    subplot(nr,nc,2)
    plot(abs(bb{ii}));grid on %<-- plot pulses for last gradient ...
    subplot(nr,nc,[3 4])
    plot(G);grid on
    hold on
    plot(Gv{ii},'--')
    hold off
end



%%



figure(1)
clf
nr=4;nc=2;
subplot(nr,nc,1)
imagesc(abs([Psol(:,:,sl)]))
axis off
title('Absolute exicitation')

subplot(nr,nc,2)
plot(tv*1e3,1e3*abs(bb{end}));
title('RF pulse amplitudes')
grid on 
xlim([tv(1) tv(end)]*1e3)
xlabel('ms')
ylabel('uT')

tag={'Gx','Gy','Gz'};
for ix=1:3
    subplot(nr,nc,[1 2]+ix*2)
    plot(t*1e3,G(:,ix),'k')
    grid on
    hold on
    plot(tv*1e3,Gv{end}(:,ix),'r--')
    xlim([0 12.1])
    legend(['Initial ' tag{ix}],['Final ' tag{ix}],'location','eastoutside')
    ylabel('mT m^{-1}')
    xlabel('ms')
end
scsz = get(0,'ScreenSize');
set(gcf,'position',[0.1*scsz(3) 0.1*scsz(4) 0.4*scsz(3) 0.6*scsz(4)])

figure(2)
clf
nr=3;nc=2;
subplot(nr,nc,1)
imagesc(imag(Psol(:,:,fix(n3/2))),[0 pi/2])
axis off
title('CPMG, axial')
subplot(nr,nc,2)
imagesc(real(Psol(:,:,fix(n3/2))),[-pi/18 pi/18])
axis off
title('non-CPMG, axial')
subplot(nr,nc,3)
imagesc(squeeze(imag(Psol(:,fix(n2/2),:))),[0 pi/2])
axis off
title('CPMG, sagittal')
subplot(nr,nc,4)
imagesc(squeeze(real(Psol(:,fix(n2/2),:))),[-pi/18 pi/18])
axis off
title('non-CPMG, sagittal')
subplot(nr,nc,5)
imagesc(squeeze(imag(Psol(fix(n1/2),:,:))),[0 pi/2])
axis off
title('CPMG, coronal')
subplot(nr,nc,6)
imagesc(squeeze(real(Psol(fix(n1/2),:,:))),[-pi/18 pi/18])
axis off
title('non-CPMG, coronal')
scsz = get(0,'ScreenSize');
set(gcf,'position',[0.5*scsz(3) 0.1*scsz(4) 0.4*scsz(3) 0.8*scsz(4)])
