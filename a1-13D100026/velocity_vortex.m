function [u,v] = velocity_vortex(x,y,x0,y0,x1,y1,Q0,Q1)
    u = Q0/(2*pi)*(y-y0)/((x-x0)^2 + (y-y0)^2) + Q1/(2*pi)*(y-y1)/((x-x1)^2 + (y-y1)^2)
    v = Q0/(2*pi)*(x-x0)/((x-x0)^2 + (y-y0)^2) + Q1/(2*pi)*(x-x1)/((x-x1)^2 + (y-y1)^2)
end
