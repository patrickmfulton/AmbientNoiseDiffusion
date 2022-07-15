function [G,F,Shift]=f_synth_ambdiff_freq_PandFlow(D,t,x,Shift)
% x is a scalar.  index is frequencies
% Returns synthetic Green's function response for diffusion in frequency domain
% with random phase applied. G is green's function for Pressure, and F is
% for dP/dz. 

i=sqrt(-1);

N=length(t);
dt=t(2)-t(1);


f=(0:N-1)./(N*dt);
w=2*pi*f;

if ~exist('Shift','var')
    MaxShift=max(t)*2;
    Shift=(rand(1)-0.5)*MaxShift;
end

%%%%GREEN'S FUNCTION
alpha=(1+i).*sqrt(w./(2*D));
if (x<0)
    G=(1./(2*D*alpha)).*exp(alpha.*x); %Shamsalsadati & Weiss, GEOPHYSICS, 2014
else
    G=(1./(2*D*alpha)).*exp(-alpha.*x); %Shamsalsadati & Weiss, GEOPHYSICS, 2014
end

G(f==0)=0;
% implementation of random phase
G=G.*exp(i*Shift.*w);


%%%%spatial derivative of Green's Function
if (x<0)
    Fo=(1./(2*D*alpha)).*exp(alpha.*x); %Shamsalsadati & Weiss, GEOPHYSICS, 2014
    F=-Fo.*alpha; %spatial derivative
else
    Fo=(1./(2*D*alpha)).*exp(-alpha.*x); %Shamsalsadati & Weiss, GEOPHYSICS, 2014
    F=-Fo.*alpha; %spatial derivative
end

F(f==0)=0;
% implementation of random phase
F=F.*exp(i*Shift.*w);
