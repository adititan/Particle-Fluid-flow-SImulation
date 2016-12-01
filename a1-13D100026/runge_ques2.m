U = 1;         
Q0 = 1;
Q1 = -1;
hx = 0.05;         
hy = 0.05;         
l=linspace(-2,2,10);
x0=-1; y0=0;
x1=1; y1=0;          
figure;

for j=1:10
x_0=-2;
y_0= l(j);
	for i=1:100
		[u,v] = velocity(x_0,y_0,x0,y0,x1,y1,U,Q0,Q1);
		k1x = hx*u;
		k1y = hy*v;
        [u0,v0] = velocity(x_0+hx,y_0+k1x,x0,y0,x1,y1,U,Q0,Q1);
		k2x = hx*u0;
		[u1,v1] = velocity(x_0+hy,y_0+k1y,x0,y0,x1,y1,U,Q0,Q1);
        k2y = hy*v1;
        x_1 = x_0 + 0.5*(k1x+k2x);
		y_1 = y_0 + 0.5*(k1y+k2y);
		arrowline([x_0,x_1],[y_0,y_1]);
		x_0 = x_1;
		y_0 = y_1;
	end
end
