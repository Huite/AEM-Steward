clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.2.4, Figure 1.30b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// Exponential integral, in equation (1.62) and from Abramowitz and Stegun (1972)  
function [out]=expint(u)
  c0=-0.57721566;c1=0.99999193;c2=-0.24991055;c3=0.05519968;c4=-0.00976004;c5=0.00107857;
  a1=8.57332897401;a2=18.0590169730;a3=8.6347608925;a4=0.2677737343;
  b1=9.5733223454;b2=25.6329561486;b3=21.0996330827;b4=3.9584969228;
  out=(-log(u)+c0+c1*u+c2*u.^2+c3*u.^3+c4*u.^4+c5*u.^5).*(abs(u)<=1)...
     +((u.^4+a1*u.^3+a2*u.^2+a3*u+a4)./((u.*exp(u)).*(u.^4+b1*u.^3+b2*u.^2+b3*u+b4))).*(abs(u)>1);
endfunction

xmin=-0.1; xmax=0.1; ymin=-0.1; ymax=0.1;
xrange =xmin+(xmax-xmin)*[0:1/400:1];   yrange= ymin+(ymax-ymin)*[0:1/400:1]';   x=(xrange'*ones(yrange'))';    y=(ones(xrange' )*yrange')';
xrangev=xmin+(xmax-xmin)*[1/20:1/10:1]; yrangev=ymin+(ymax-ymin)*[1/20:1/10:1]'; xv=(xrangev'*ones(yrangev'))'; yv=(ones(xrangev')*yrangev')';
r=sqrt(x.^2+y.^2); r0=1e-6; r=r+(r0-r).*(r<r0); rv=sqrt(xv.^2+yv.^2); thetav=atan(yv,xv);
diffusivity=100; t=10^(-4); Q=1;  Phi0=0;
u=r.^2/(4*diffusivity*t);
Phi=Phi0-Q/(4*%pi)*expint(u);
uv=rv.^2/(4*diffusivity*t);
qx=-Q/(2*%pi)*exp(-uv)./rv.*cos(thetav); qx=qx.*(rv>0.02);
qy=-Q/(2*%pi)*exp(-uv)./rv.*sin(thetav); qy=qy.*(rv>0.02);
draw_flow_conduction(xrange,yrange,Phi,[-1:.05:0],[],[],xrangev,yrangev,qx,qy);
xtitle(title_author,"x","y");
draw_points(0,0);
