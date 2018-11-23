function [F] = getForce(c00, c20, E, nu)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
mu=E/(2*(1+nu));
%N=(4*nu*(2-3*nu)+1-5*(1-2*nu)+12*(2-3*nu)*(1-2*nu))/(4*(2-3*nu)*(1-2*nu));
%F=(mu*c20*sqrt(5)*c00*N)/4;
N=(2*nu/(1-2*nu)-(1+1/(1-2*nu))/(2*(2-3*nu))+1);
F=(mu*c20*sqrt(5)*c00*N)/2;
end

