Q0=2*pi;
Q1=2*pi;          
xmin = -1;       
xmax = 2;        
dx = .11;         
ymin = -1;        
ymax = 1;         
dy = .12;         

gcontf = 25;               
gcont = 25;                                      
gstream = 1;                

zstart0 = 0 + i*0;
zstart1 = 1 + i*0;

[x,y] = meshgrid ([xmin:dx:xmax],[ymin:dy:ymax]);  
z = x+i*y; 
Phi =  -1*i/(pi+pi)*[Q0*log(z-zstart0) + Q1*log(z-zstart1)];   
      

if gcontf
    colormap(winter);  
    contourf (x,y,real(Phi),gcontf);
end
hold on;
if gcont contour (x,y,imag(Phi),gcont,'w'); end
