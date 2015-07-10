% Example 2D spiral trajectory phase relaxed 90° pulse.
% Shaihan Malik 9th July 2015. 

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

%% Load in B0 and B1 maps, plus mask

% B0 is in units of uT. B1 field ('tx') is relative sensitivity
load b0_b1_maps_2d
[n1 n2 Nc]=size(tx);           

% Make list of coordinates, and index of those within mask
xlist = double(cat(2,X(:),Y(:),Z(:)));
idx = find(m);
Nidx=length(idx);


%% target
P = zeros([n1 n2]);

% Load in heart shape
load heart

%%% smooth heart
h = fspecial('gaussian',[20 20],4);
heart = imfilter(heart,h);

%%% Interpolate down to off-centre location
hsize = 55e-3; %% size of heart
[xi,yi] = meshgrid(linspace(0,hsize,size(heart,2)),linspace(0,hsize,size(heart,1)));
xoff = -10e-3;
yoff = 10.e-3;

P = interp2(xi+xoff,yi+yoff,heart,X,Y);
P(isnan(P))=0;
%%% apodize
P = GaussConv3D(P,2.5e-3,m,X,Y,Z);

%%% Now add pi/2 phase and scale to flip
P=P*flip*1i;

figure;
imagesc(abs(10*P+m));




%% Speed-limited k-space trajectory
usfactor=3.5;
%%% DEFINE TARGET TRAJECTORY - 2D spiral
fovx=(max(X(:))-min(X(:)));
dx=3.4e-3; %<- design resolution 3.4mm
kmax = pi/dx;
dk=usfactor*2*pi/fovx;
nturn=ceil(kmax/dk);
M=1000;
% parametric spiral generation
tau=linspace(0,1,M);
tf = 1./(1+exp(30*(tau-0.9)));
tf = cumsum(tf);
tf = tf/max(tf);
kr=kmax*(1-tf);
th = 2*pi*nturn*tf;
K = [kr(:).*cos(th(:)) kr(:).*sin(th(:)) 0*th(:)];%<-- X-Y plane
%%% explicitly set last value to zero

%%% explicitly set last value to zero
K = cat(1,K,zeros([1 3]));
M=M+1;

%%% First gradient trajectory - this only uses gradient constraints (not RF)
% This is Lustig's time optimal gradient code (DOI:10.1109/TMI.2008.922699)
Kcm = K / (2*pi*100);% k-space in 1/cm
Gmax = 30 / 10; % G/cm
Smax = 100 / 10; % G/cm/ms

%%% Make simple initial gradient to get k-space velocity weighting
[~,G] = gradient(Kcm/4.257e3,dt);

%%% Get velocity limited k-space trajectory by using VERSE guided optimal
%%% k-space (doi: 10.1002/mrm.21950) with fictional RF pulse whose peak amplitude 
%%% ramps up towards 25uT at the end. The VERSE guided design is used to 
%%% slow the trajectory at the end to keep this peak at 10uT - i.e. the 
%%% trajectory linearly decelerates towards the end.
b1max = 10; % This is in uT
blimit = ones([M 1])*1e-3;
blimit(fix(M*0.6):end)=linspace(1e-3,25e-3,fix(0.4*M)+2);

opt = reVERSE_init;
opt.b = blimit*10;%<-- mT->Gauss
opt.Gmod = sum(G.^2,2).^0.5;
opt.bmax = b1max*1e-2;%<-- uT->Gauss
[Cv,time,g,s,k, phi, sta, stb,p_of_t,oo] = minTimeGradient_VERSE(Kcm,0,0,Gmax,Smax,dt*1e3,[],[],opt);

fprintf(1,'Duration = %1.3f ms\n',time);

%%% check actual TRANSMIT kspace
G = 10*g; %<- convert back to mT/m
M=length(G);

%%% simulate the effect of Eddy currents, use measured GIRF from KCL 3T
%%% 8-ch PTX system
girf = load('GIRF_3T_London_20140729');
gcor = @(x)(gradient_distort_GIRF(x,girf.ff,girf.Hw,dt,10));

[Gedd,kedd] = gcor(G);

%%% add time if using multi-frequency
t = (1:M)*dt;
kt = t(:)-t(end);
klist=kedd;



%% Work out L-curves using mcgls. Remember that the error doesn't really 
%%% change much when VERSE is applied, so can get the starting point here
%%% MCGLS is available from: http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=3&SC_ID=10&MP_ID=126
%%% and is published here: doi:10.1016/j.apnum.2003.11.010
%%% If MCGLS is not available then use another linear matrix inversion
%%% routine here.

L = [0.1 1 8:32:96 128:64:512 640:128:1024]; %<-- Array of Lambda values
Wreal = [0 0.1 0.2 0.5 1];                   %<-- Array of weightings for non-CPMG

nl = length(L);
nw=length(Wreal);

%%% Define system matrix once
fprintf(1,'Define STA system matrix\n');
A = STA_system_matrix(xlist,klist,t,tx,b0,m,'loopcalc');
% Apply fat weighting
wA = A;
fprintf(1,'Split RE and IM\n');
%%% Now split real and imaginary parts
Aw = [imag(wA) real(wA);real(wA) -imag(wA)];
Pw = [imag(P(idx));real(P(idx))];

bb={};
Psol={};
err=[];pwr_av=[];pwr_pk=[];
%%% Loop over Wreal weightings
for widx=1:nw
    % add a weighting in here for the real/imag channels
    WW = (Wreal(widx)+1i)*double(m);
    
    ww = [imag(WW(idx));real(WW(idx))];
    Aw_tmp = Aw .* repmat(ww,[1 size(Aw,2)]);
    Pw_tmp = ww.*Pw;
    
    fprintf(1,'W value %d: Solve LLS problem with mCGLS\n',widx);
    opt=struct('max_it',50);
    tic;
    [b,brelres] = mcgls(Aw_tmp,Pw_tmp,L,opt);
    toc
    fprintf(1,'... done. Relative residual %1.3f, max B1 %1.2f uT\n',brelres(end,1)/norm(Pw(:)),1000*max(abs(b(:))));
    
    
    %%% now look at excitation pattern
    b = complex(b(1:M*Nc,:),b(M*Nc+1:M*Nc*2,:));
    clear Aw_tmp Pw_tmp
    
    % put into real space
    Psol{widx} = zeros([prod(size(P)) nl]);
    Psol{widx}(idx,:) = A * b(:,:);
    Psol{widx}=reshape(Psol{widx},[n1 n2 nl]);
    
    % look at errors in real and imaginary channels
    for ii=1:nl
        Ptmp = Psol{widx}(:,:,ii);
        err(ii,widx,1)=norm(real(P(:))-real(Ptmp(:)))/norm(P(:));
        err(ii,widx,2)=norm(imag(P(:))-imag(Ptmp(:)))/norm(P(:));
        %err(ii,jj,3)=norm(PP(:)-Pw(:))+norm(L(ii)*bb{ii,jj}(:));
        pwr_av(ii,widx)=1e3*norm(b(:,ii)).^2;
        pwr_pk(ii,widx)=max(abs(b(:,ii)));
    end
    
    b = reshape(b,[M Nc nl]);
    bb{widx}=b;
end    

%%% Find solution indices closest to peak power limit
solix = ones([nw 1]);
max_b1 = 13; % uT
for ii=1:nw
    tmp=min(find(pwr_pk(:,ii)<max_b1*1e-3));

    if ~isempty(tmp)
        solix(ii)=tmp;
    else
        % no solution so pick highest lambda
        solix(ii)=nl;
    end
end




%% Display the results

figure(1)
nr=2;nc=2;
cmap=colormap(hsv);
pp=zeros([nw 2]);
titls={'Imag (CPMG)','Real (non-CPMG)'};
xl={[0 0.55],[0 1]};
for ii=1:2
    subplot(nr,nc,ii)
    pp(:,ii)=plot(err(:,:,3-ii),pwr_av);
    grid on
    title(titls{ii},'fontsize',12)
    ylabel('RF power')
    xlabel('nRMSE')
    hold on
    for jj=1:nw
        plot(err(solix(jj),jj,3-ii),pwr_av(solix(jj),jj),'or')
    end
    
    subplot(nr,nc,ii+nc)
    pp(:,ii+2)=plot(err(:,:,3-ii),pwr_pk);
    grid on
    title(titls{ii},'fontsize',12)
    ylabel('max B1')
    xlabel('nRMSE')
    hold on
    for jj=1:nw
        plot(err(solix(jj),jj,3-ii),pwr_pk(solix(jj),jj),'or')
    end
end
scsz = get(0,'ScreenSize');
set(gcf,'position',[0.5*scsz(3) 0.2*scsz(4) 0.4*scsz(3) 0.6*scsz(4)])


leg={};for ii=1:nw,leg{ii}=sprintf('W=%1.1f',Wreal(ii));end
legend(leg,'location','northeast')


%%%

figure(2);
nr=3;nc=nw;
for ii=1:nw
    subplot(nr,nc,ii)
    imagesc(imag(Psol{ii}(:,:,solix(ii))),[-0.1 1.1]*flip)
    title(sprintf('CPMG W=%1.2f L=%d',Wreal(ii),fix(L(solix(ii)))))
    axis off
    
    subplot(nr,nc,ii+nc)
    imagesc(real(Psol{ii}(:,:,solix(ii))),[-0.5 0.5]*flip)
    title(sprintf('nonCPMG W=%1.2f L=%d',Wreal(ii),fix(L(solix(ii)))))
    axis off
    
    subplot(nr,nc,ii+nc*2)
    imagesc(abs(Psol{ii}(:,:,solix(ii))),[0 1.1]*flip)
    title(sprintf('abs. W=%1.2f L=%d',Wreal(ii),fix(L(solix(ii)))))
    axis off
    
end
set(gcf,'position',[0.1*scsz(3) 0.15*scsz(4) 0.6*scsz(3) 0.5*scsz(4)])


