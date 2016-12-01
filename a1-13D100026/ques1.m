U = 1;         
Q = 1;           
xmin = -2.5;       
xmax = 2.5;        
dx = .11;         
ymin = -2.5;        
ymax = 2.5;         
dy = .12;         

gcontf = 25;               
gcont = 25;                                      
gstream = 1;                

zstart = xmin + i*(ymin+(ymax-ymin)*[0.1:0.1:1]);
zstart = [zstart xmin-i*0.001];

z0=1+i*0;
[x,y] = meshgrid ([xmin:dx:xmax],[ymin:dy:ymax]);  
z = x+i*y; 
Phi = U*z + (Q/(pi+pi))*[log(z+z0)-log(z-z0)];   
      

if gcontf
    colormap(winter);  
    contourf (x,y,real(Phi),gcontf);
end
hold on;
if gcont contour (x,y,imag(Phi),gcont,'w'); end
