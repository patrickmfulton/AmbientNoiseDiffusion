function [dGm]=f_DerDiffGreen(t,L,D,dim)
% [dGm]=f_DerDiffGreen(t,L,D,dim)
% Gives the spatial derivative (dz) of the (1D, 2D, or 3D) Green's function 
% response for diffusion, dGm, as a function of time, t at distance L (i.e. 
% dz from source) away for diffusivity D, assuming dx and dy = 0. 
% Assumes delta function source, but output can be scaled by a source 
% amplitude.
% 
% IMPORTANT!!!
% Note: Inputs must have comparable units [e.g. t[day], L[m], D[m^2/day]
%
% INPUTS:
% t (vector) [time]
% L (constant or vector) [length]
% D (constant or vector) [length^2 / time]
% dim (integer) [1=1D 2=2D 3=3D] Dimensions of diffusion
%
% Can be run for multiple L or D as a vector, but not both.
%
% Patrick Fulton
% pfulton@cornell.edu
% March 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<4 %if no choice for 1D, 2D, or 3D, assume 1D
    dim=1;
end


if size(t,2)>size(t,1)
t=t';
end
if size(L,1)>size(t,2)
L=L';
end
if size(D,1)>size(t,2)
D=D';
end

if dim == 1

   %%%%%% 1D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [(2L)/(4Dt)] * [1/(2(pi D t]^(1/2)) * [exp(-(L^2)/(4 D t)]  
   dGm=(2.*L./(4.*D.*t)).*1./(2*(pi.*D.*t).^(1/2)).*exp(-(L).^2./(4*D.*t));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif dim == 2

   %%%%%% 2D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [(2L)/(4Dt)] * [1/(4(pi D t]^(2/2)) * [exp(-(L^2)/(4 D t)]  
   dGm=(2.*L./(4.*D.*t)).*1./(4*(pi.*D.*t).^(2/2)).*exp(-(L).^2./(4*D.*t));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elseif dim == 3

   %%%%%% 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [(2L)/(4Dt)] * [1/(8(pi D t]^(3/2)) * [exp(-(L^2)/(4 D t)]  
   dGm=(2.*L./(4.*D.*t)).*1./(8*(pi.*D.*t).^(3/2)).*exp(-(L).^2./(4*D.*t));
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



    
    
    
    
    
    
    
    