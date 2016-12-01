Q0=2*pi;
Q1=2*pi;
d=1; %distance
r=0.5*d;
x0=0; y0=0;
x1=1; y1=0;
xm= 0.5*(x0+x1);
ym= 0.5*(y0+y1);
h=0.01;
n =400;
figure;
[u,v] = vel1_dueto2(x1,y1,x0,y0,Q0);
x_0 =x1;
y_0 =y1;
for j=1:n
        [u1,v1] = vel_circle(xm,ym,x_0,y_0,v);
        x_1 = x_0 + h * u1;
		y_1 = y_0 + h * v1;
        arrowline([x_0,x_1],[y_0,y_1]);
		x_0 = x_1;
		y_0 = y_1;
end
time=h*n;
time
theta=(v*time)/r;
x=r*cos(theta);
y=r*sin(theta);
dist_error=sqrt((x_0-x)^2+(y_0-y)^2)