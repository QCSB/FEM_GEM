function [ mu ] = myfunGenome_growth_warm_xy_phase2(c,c2,x,y,t)

global Vmax
global Km
global growth
global Yield
global V
if y>0.08
    mu=0;
else
mu=growth;
end
end