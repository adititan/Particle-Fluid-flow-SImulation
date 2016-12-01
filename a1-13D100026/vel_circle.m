function [u1,v1] = vel_circle(xm,ym,x1,y1,v)
    r=sqrt((xm-x1)^2+(ym-y1)^2);
    u1=-v*(y1-ym)/r;
    v1=v*(x1-xm)/r;

end
