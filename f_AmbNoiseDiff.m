function [tw,G,dlag]=f_AmbNoiseDiff(t,X1,X2,win,dw,scaleopt)
% [tw,G,dlag]=f_AmbNoiseDiff(t,X1,X2,win,dw,scaleopt)
%
%Perform ambient noise diffusion analysis with time series from two
% observation points and extract response curve, G, as a function of lag.

% win is window size
% dw = # of samples to window advance each step;
% OUTPUT:
% tw: max time of window (causal), 
% G: -2diff(xcorr(X1,X2) (the response curve; sometimes equal to Greens fn)
% dlag: (lag)
%
%
%Patrick Fulton | pfulton@cornell.edu | April 2017 | updated July 2022

if nargin~=6
    scaleopt='none'; %unbiased?
end

%% Make the correct dimensions
if size(X1,2)>size(X1,1)
   X1=X1';
end

if size(X2,2)>size(X2,1)
   X2=X2';
end

if length(X1)~=length(X2)
    error('Error: lengths of d1 and d2 do not match')
end
%%

dt=mode(diff(t));
N=length(t);

if win == N
    nwins =1;
else
    nwins=ceil((N-win)/dw); %round up
end


% dlag=diff(lag)/2+lag(1:end-1);
C=zeros(1,nwins);
lag=[-win+1:1:win-1]*dt;
CC=zeros(length(lag),nwins);    

parfor jj=1:nwins
    winstart = 1+(jj-1)*dw; 
    winend = winstart + win-1;
    
    X1sub=detrend(X1(winstart:winend));%,'constant');
    X2sub=detrend(X2(winstart:winend));%,'constant');
    
    
    corf=xcorr(X1sub,X2sub,scaleopt);

    
    CC(:,jj)=corf;

    
    tw(jj)=max(t(winstart:winend));
   
end

G=-2*diff(CC)/dt;

dlag=lag(1:end-1)+dt/2;

