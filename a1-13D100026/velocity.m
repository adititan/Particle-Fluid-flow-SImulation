function[u,v] = velocity(x,y,x0,y0,x1,y1,U,Q0,Q1)
r0=sqrt((x-x0)^2 + (y-y0)^2);
r1=sqrt((x-x1)^2 + (y-y1)^2);
u = U + Q0/(2*pi*r0^2)*(x-x0) + Q1/(2*pi*r1^2)*(x-x1);
v= Q0/(2*pi*r0^2)*(y-y0) + Q1/(2*pi*r1^2)*(y-y1);

end