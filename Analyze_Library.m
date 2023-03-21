clear, close all; clc

% processes synthetic library
for ii = 1:41
    for nset=1:2,
    load sim_20210514

    N=length(sim.t(ii,:));
    N2=floor(N/2);
    if nset==1
        Iseg=[1:floor(N*0.6)]; % Library
    else
        Iseg=[floor(N*0.6)+1:N]; % Test set
    end
    t=sim.t(ii,Iseg);
    data1=sim.Tbres(ii,Iseg);
    data2=sim.Tares(ii,Iseg);

    
dt=mode(diff(t));

% Do Ambient Diffusion analysis (i.e. d/dt of xcorr)
win=length(data1)/1;
[tw,Gd,dlag]=f_AmbNoiseDiff(t,data1,data2,win,win,'biased');

% Gd=Gd(round(length(Gd)/2)+1-Nav:end,:); %just positive lags (and Nav for averaging)
% dlag2=(1-Nav:length(Gd)-Nav)*dt;  %convert lags units to time
Gd=Gd(round(length(Gd)/2)+1:end,:); %just positive lags 
dlag2=(1:length(Gd))*dt;  %convert lags units to time

Gdn=Gd/max((Gd(dlag2<864000))); %normalize by max within short lag
Gdn=Gdn/Gdn(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
ax(2)=subplot(3,1,2);

plot(t/86400,[data1; data2],'linewidth',2)
xlabel('Time (days)')
ylabel('Synth. Data')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
hold on

ax(3)=subplot(3,1,3);

plot(dlag2/86400,Gdn,'-','linewidth',2), hold on

hold on
xlim([0 5])
xlabel('Lag (days)')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)

if(nset==1)
    GG(:,ii)=Gd;
    dlag2_Library=dlag2;
else
    GG_test(:,ii)=Gd;
    dlag2_test=dlag2;
end


    end
end
%%
%close all;

figure(4);
ax(1)=subplot(2,1,1);

ni=12;  %12
nj=33;  %33

nn=nj-ni+1;

set(gca, 'ColorOrder', colormap(jet(nn)), 'NextPlot', 'replacechildren');
II=ni:nj;
offsetArray=II'*ones(1,length(sim.Tares))/7e4;
plot(sim.t(1,:)/86400,[sim.Tares(ni:nj,:)]+offsetArray,'linewidth',1)
hold on
plot(sim.t(1,:)/86400,[sim.Tbres(ni:nj,:)]+offsetArray-5e-4,'linewidth',1)
errorbar([800 800],[-4.5e-4 -3.5e-4],[0 0],[0 0],[5 5] ,[5 5],'k','LineWidth',3) % scalebar
text(810,-4e-4,'10^{-4} ^o C')
hold off
yticks('');


xlabel('Time (days)')
caxis([log10(sim.D(ni)) log10(sim.D(nj))])
ylim([-5 5]*1e-4);

set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
hold on
colorbar 
box on
xlim([0 900])


figure(3)
ax(2)=subplot(2,1,2);
set(gca, 'ColorOrder', colormap(jet(nn)), 'NextPlot', 'replacechildren');
dlag2=dlag2_Library;
plot(dlag2/86400,GG(:,ni:nj),'-','linewidth',2), hold on  %12:34

  ylim([-7 1]*1e-16)
  caxis([log10(sim.D(ni)) log10(sim.D(nj))])
xlim([0 1])
xlabel('Lag (days)')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
    box on
colorbar 


%%
%matching code -finds best fit from library for each synthetic from 1 day
%of lags
Nav=3; % smoothing average
L=1.5; % sensor spacing

ndays_max=4;
Itime=(dlag2<86400*ndays_max); % put some extra point for smoothing
Library=GG(Itime,ni:nj);
Test=GG_test((dlag2_test<86400*ndays_max),ni:nj);
Nlibrary=nj-ni+1;
% smooth the library
Library_smooth=[];
Test_smooth=[];
for i=1:Nlibrary,
    Library_smooth(:,i)=conv(Library(:,i),ones(1,Nav)/Nav,'same');
    Library_smooth(:,i)=Library_smooth(:,i)./range((Library_smooth(:,i)));

    Test_smooth(:,i)=conv(Test(:,i),ones(1,Nav)/Nav,'same');
    Test_smooth(:,i)=Test_smooth(:,i)./range((Test_smooth(:,i)));

end
% now cut down to 1 day
Itime=(dlag2<86400*1)&(dlag2>00);
Library_smooth=Library_smooth(Itime,:);
Test_smooth=Test_smooth(Itime,:);

D=sim.D(ni:nj);

Ntest=nj-ni+1;
for itest=1:Ntest;
test=Test_smooth(:,itest);

RMS=sqrt(sum((Library_smooth-test).^2));
[Res(itest),Ibest]=min(RMS);
% 
% [D(Ibest) D(itest)]
% [Ibest itest]
%Fit_Diff(itest)=Ibest-itest;
Fit_D(itest)=D(Ibest);

figure(7)
plot(test);
hold on
plot(Library_smooth(:,Ibest))
hold off
pause
end

figure(8)
loglog(D,Fit_D,'r*');
xlabel('Modeled D (m^2/s)')
ylabel('Best Fit D (m^2/s)')
set(gca,'FontSize',18)

%PercentFit=sum(abs(Fit_Diff)<=2)/length(Fit_Diff)
