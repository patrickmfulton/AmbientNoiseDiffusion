% MakeFigure_DIffusionCurves
% Compute analytical solution for the impulse response for 1,2,or 3D 
% diffusion to demonstrate the effect of diffusivty on the time-series
% response observed at a given distance away from the source.
%
% Constructs figure that forms the basis of Figure 1 in manuscript by
% Fulton & Brodsky: "Determining Hydraulic Diffusivity from Ambient Noise 
% Signals of Advection in Temperature Data"
%
% mks units used throughout for calculations. 
% Time plotted in days for illustrative purposes.
%
% Utilizes custom functions: 
% f_DiffGreen
% save2pdf (optional)
%
% Patrick M. Fulton
% pfulton@cornell.edu
% July 2022
% https://github.com/patrickmfulton/AmbientNoiseDiffusion

%%
clear, close all; clc

t=[0:.005:20]*86400; %time vector (seconds)
n=1; % dimensions for diffusion 1, 2, or 3.
z=[0 1.5]; %distance between source and observation point
D=10.^[-4:-.5:-5]; %Diffusivity [m2/s]

for zz=1:length(z) 
   
for dd=1:3
    
    [Gm_tmp]=f_DiffGreen(t,z(zz),D(dd),n); %Analytical soln.

    A(:,dd,zz)=Gm_tmp; %Save response curves for various D or z values
    
end
    
end
    
%% Plot results

maxA=max(A(:,:,1));

figure(1)
hold on
plot(t/86400,A(:,1,2)./maxA(1),'-k','linewidth',3) % z=1.5 (or non-zero) scenario
plot(t/86400,A(:,2,2)./maxA(2),'-b','linewidth',3) % z=1.5 (or non-zero) scenario
plot(t/86400,A(:,3,2)./maxA(3),'-c','linewidth',3) % z=1.5 (or non-zero) scenario
legend(['D=1E',num2str(log10(D(1)))],['D=1E',num2str(log10(D(2)))],['D=1E',num2str(log10(D(3)))])


xlim([0 2])
xlabel('Time (days)')
ylabel('Amplitude (m)')
set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
    hold on
    box on
    axis square


xlim([0 2])
xlabel('Time (days)')
ylabel('Amplitude (m)')
title('Effect of Diffusivity on Impulse Response (z=1.5 m)')

set(gca,'fontsize',16)
    set(findall(gcf,'type','text'),'fontsize',16)
    hold on
    box on
    axis square
    
    %% Save figure
%   save2pdf('Figure1_Diffusivity_ImpulseCurves_test')
% %https://www.mathworks.com/matlabcentral/fileexchange/16179-save2pdf

