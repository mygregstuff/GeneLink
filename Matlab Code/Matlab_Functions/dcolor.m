function [out]=dcolor(n,varargin)
out=distinguishable_colors(n);
set(0,'DefaultAxesColorOrder',out)
end