%% AmbientNoiseDiffusion_synthetics.m
%
% Demonstration of ambient noise diffusion analysis as illustrated in
% Figure 2 of manuscript by Fulton & Brodsky: "Determining Hydraulic 
% Diffusivity from Ambient Noise Signals of Advection in Temperature Data"
% which builds upon the work of Sneider, 2006; Fan and Sneider, 2009; 
% and Shamsalsadati and Weiss, 2014.
%
% Objectives:
% 1) Create synthetic ambient noise time-series of pressure and pressure 
% gradient at 2 observation points resulting from distributed pressure 
% sources.
% 2) perform ambient noise diffusion analysis to reconstruct Green's
% function response
% 3) compare empirical responses to analytical solutions 
% mks units used throughout for calculations. 
% Time or lag vectors plotted in days for illustrative purposes.
%
% Utilizes the following custom functions:
% f_AmbNoiseDiff - uses parfor parallelization; 
%       edit 'parfor' to 'for' to unparallelize if necessary.
% f_DiffGreen
% f_DerDiffGreen
% f_synth_ambdiff_freq_PandFlow
% save2pdf (optional)
%
%
% Patrick M. Fulton
% pfulton@cornell.edu
% July 2022
% https://github.com/patrickmfulton/AmbientNoiseDiffusion


clear, close all; clc

%%
savepdf=0; %save pdf?
R=1;% Number of runs to average over.

%%
%%MODEL PARAMETERS    

dD=.1;
% Dall=10.^(-10:dD:0); %initial range of D estimates
Dall=10.^-(5); % option to perform on single D estimate

L=1.5;  % distance between two receivers [m];

dt=600;
tmax=86400*900;

t=[dt:dt:tmax]; %time vector [s]

%%SOURCE DISTRIBUTION
Ws=34; %width of source distribution
rhos=1.147; % Source density

%location of two receivers (i.e. observation points / sensors)
x_R1=L/2;
x_R2=-L/2; % symetric around zero

%location of source 
dx_source=1./rhos;
x_source1=dx_source/2:dx_source:Ws/2;
x_source2=(-dx_source/2:-dx_source:-Ws/2);

x_source=[x_source2 x_source1];

% Plot spatial distribution of sources and observation points
figure(1)
ax(1)=subplot(1,5,1);
plot(1,[x_R1 x_R2],'k>')
hold on
plot(1,[x_source2 x_source1],'r*')
ylim([-1 1]*20)
axis ij
xlabel('z-position (m)')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
    
%%
% SOLUTION OF PULSE SOURCE



N=length(t);
f=(0:N-1)/(N*dt); %frequency
w=2*pi*f; %angular frequency
GG=zeros(length(t)-1,R); %initialize array for Green's function responses for Pressure
FF=zeros(length(t)-1,R); %initialize array for Green's function responses for Pressure gradient



for dd=1:length(Dall)
    D=Dall(dd);
for ii=1:R % RUN SEVERAL TIMES    
    
    
%%Analytical Green function between two receiver
G_th=f_DiffGreen(t,L,D,1)'; %theoretical soln for Pressure
F_th=f_DerDiffGreen(t,L,D,1)'; %theoretical soln for pressure gradient
 


%% Create synthetic signals for analysis
% Diffusion signal at x1 from source ix in freq. domain with random phase
% and Diffusion signal at x2 from source ix in freq. domain...
%     maintaining phase shifts for both receivers, x1 and x2.
%  G represents Presure response; F represents pressure gradient response.

for ix=1:length(x_source)

 x1=(x_R1-x_source(ix));
 x2=(x_R2-x_source(ix));


[G1,F1,Shift]=f_synth_ambdiff_freq_PandFlow(D,t,x1); 
[G2,F2]=f_synth_ambdiff_freq_PandFlow(D,t,x2,Shift); 

G1A(ix,:)=G1; %matrix of all the source responses for Pressure
G2A(ix,:)=G2;

F1A(ix,:)=F1; %matrix of all the source responses for Pressure gradient
F2A(ix,:)=F2;

end
%% 
%% Construct synthetic time-series by compiling source responses 
%  and transforming to time-domain.
dw=w(2)-w(1);

G1=sum(G1A,1); %synthetic series  of all sources
G2=sum(G2A,1);
data1=real(ifft(G1)*dw); %synthetic time series for R1
data2=real(ifft(G2)*dw); % " for R2

F1=sum(F1A,1); %synthetic series of all sources
F2=sum(F2A,1);
ddata1=real(ifft(F1)*dw); %synthetic time series for R1
ddata2=real(ifft(F2)*dw); % " for R2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do Ambient Diffusion analysis (i.e. -2*d/dt of xcorr)
win=length(data1)/1; %time window to evaluate over = full time series. 
[tw,Gd,dlag]=f_AmbNoiseDiff(t,data1,data2,win,win,'biased');

Gd=Gd(round(length(Gd)/2)+1:end,:); %just positive lags
dlag2=(1:length(Gd))*dt;  %convert lags units to time

Gdn=Gd/max(Gd(dlag2<864000)); %normalize by max within short lag
GG(:,ii)=Gd; %store empirical response from separate runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
win=length(data1)/1;
[tw,Fd,dlag]=f_AmbNoiseDiff(t,ddata1,ddata2,win,win,'biased');

Fd=Fd(round(length(Fd)/2)+1:end,:); %just positive lags
dlag2=(1:length(Gd))*dt;  %convert lags units to time

Fdn=Fd/max(Fd(dlag2<864000)); %normalize by max within short lag
FF(:,ii)=Fd;  %store empirical response from separate runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot synthetic time series and empirical & analytical response functions
figure(2)
ax(1)=subplot(3,2,1);
plot(t/86400,[data1; data2],'linewidth',2)
xlabel('Time (days)')
ylabel('Synth. Data')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
    
ax(2)=subplot(3,2,3);
plot(dlag2/86400,Gdn,'k-','linewidth',2), hold on

hold on

    plot(t/86400,G_th./max(G_th),'r--','linewidth',2)

xlim([0 5])
xlabel('Lag (days)')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)

ax(4)=subplot(3,2,5);
plot(dlag2/86400,Fdn,'k-','linewidth',2), hold on
hold on
plot(t/86400,F_th./max(F_th),'r--','linewidth',2)

xlim([0 5])
ylim([-1 1])
xlabel('Lag (days)')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
%%



end
end

if savepdf==1
save2pdf(['Figure2ABC_',datestr(now,'yyyymmddhhMM')])  
% https://www.mathworks.com/matlabcentral/fileexchange/16179-save2pdf 
end



