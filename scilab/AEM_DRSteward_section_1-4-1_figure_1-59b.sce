clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.4.1, Figure 1.59b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

function [sigmax,sigmay,tauxy]=cantilever_stress(x,y,Lx,Ly,F,E,nu)
  I=2/3*Ly^3;
  sigmax=-F/I*(x-Lx).*y;
  sigmay=zeros(x);
  tauxy =F/(2*I)*(y.^2-Ly^2);
endfunction
function [ux,uy]=cantilever_displacement(x,y,Lx,Ly,F,E,nu)
  I=2/3*Ly^3;
  ux=-F/(E*I)*(x.^2/2-Lx*x).*y+F*(2+nu)/(6*E*I)*y.^3-F*(1+nu)/(E*I)*Ly^2*y;
  uy=F*nu/(2*E*I)*(x-Lx).*y.^2+F/(2*E*I)*(x.^3/3-Lx*x.^2);
endfunction
  
Lx=0.5; Ly=0.05; F=1*10^6; E=200*10^9; nu=0.3;
xmin=0; xmax=Lx; ymin=-Ly; ymax=Ly;
xrange =xmin+(xmax-xmin)*[0:1/400:1];  yrange= ymin+(ymax-ymin)*[0:1/200:1]';  x=(xrange'*ones(yrange'))';    y=(ones(xrange' )*yrange')';
xrangev=xmin+(xmax-xmin)*[1/40:1/20:1]; yrangev=ymin+(ymax-ymin)*[1/10:1/5:1]'; xv=(xrangev'*ones(yrangev'))'; yv=(ones(xrangev')*yrangev')';

[sigmax,sigmay,tauxy]=cantilever_stress(x,y,Lx,Ly,F,E,nu);
[vec.sigmax,vec.sigmay,vec.tauxy]=cantilever_stress(xv,yv,Lx,Ly,F,E,nu);
[ux,uy]=cantilever_displacement(xv,yv,Lx,Ly,F,E,nu); vec.ux=ux; vec.uy=uy;
  
// visualize stress
draw_stress(xrange,yrange,sigmax,sigmay,tauxy,[-2:.2:2]*10^8,xrangev,yrangev,vec,%T);
// visualize displacement
//draw_stress(xrange,yrange,sigmax,sigmay,tauxy,[0:.1:2]*10^8,xrangev,yrangev,vec,%F);

fig=gcf(); fig.figure_size(2)=370;

xtitle(title_author,"x","y");
