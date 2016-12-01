function [u,v] = vel1_dueto2(x1,y1,x0,y0,Q0)
    u = Q0/(2*pi)*(y1-y0)/((x1-x0)^2 + (y1-y0)^2);
    v = -Q0/(2*pi)*(x1-x0)/((x1-x0)^2 + (y1-y0)^2);
end
