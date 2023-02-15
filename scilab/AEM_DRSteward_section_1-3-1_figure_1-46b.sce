clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.3.1, Figure 1.46b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

T=4; h=10; A0=0.4; g=9.80665; omega=2*%pi/T;

// compute wave number
k=omega^2/g; kprev=0; 
while (abs(k-kprev)/k)>1e-12 do
  kprev=k;  k=kprev+(omega^2-kprev*g*tanh(kprev*h))/(kprev*g*h/((cosh(kprev*h))^2)+g*tanh(kprev*h));
end
L=2*%pi/k;

xmin=0; xmax=L; zmin=-h; zmax=A0;
xrange =xmin+(xmax-xmin)*[0:1/500:1];   zrange= zmin+(zmax-zmin)*[0:1/200:1]';  x=(xrange'*ones(zrange'))';    z=(ones(xrange' )*zrange')';
xrangev=xmin+(xmax-xmin)*[1/20:1/10:1]; zrangev=zmin+(zmax-zmin)*[1/10:1/5:1]'; xv=(xrangev'*ones(zrangev'))'; zv=(ones(xrangev')*zrangev')';

// free surface elevation
t=0.0*T;
surface =A0*cos(k*x -omega*t);
surfacev=A0*cos(k*xv-omega*t);

// Potential and vector components
Phi=-g*A0/(omega*cosh(k*h))*cosh(k*(z+h)).*sin(k*x-omega*t);   Phi=Phi.*(z<surface)-10*(z>=surface);
vx=g*k*A0/(omega*cosh(k*h))*cosh(k*(zv+h)).*cos(k*xv-omega*t); vx=vx.*(zv<surfacev);
vz=g*k*A0/(omega*cosh(k*h))*sinh(k*(zv+h)).*sin(k*xv-omega*t); vz=vz.*(zv<surfacev);

draw_flow_conduction(xrange,zrange,Phi,[-3:.25:3],[],[],xrangev,zrangev,vx,vz);
xtitle(title_author,"x","z");
