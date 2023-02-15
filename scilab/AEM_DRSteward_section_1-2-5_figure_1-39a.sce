clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.2.5, Figure 1.39a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

zc=1+%i; r0=2; Ex0=1; 
xmin=real(zc)-1.5*r0; xmax=real(zc)+1.5*r0; ymin=imag(zc)-1.5*r0; ymax=imag(zc)+1.5*r0;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/24:1/12:1]; yv=ymin+(ymax-ymin)*[1/24:1/12:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';
omega=-conj(Ex0)*z+Ex0*r0^2*(z-zc).^(-1); omega=omega.*bool2s(abs(z -zc)>r0)+(-conj(Ex0)*real(zc))*bool2s(abs(z -zc)<r0);
Phi=real(omega); Psi=imag(omega);
v=Ex0+conj(Ex0)*r0^2*(conj(zv-zc)).^(-2); v    =v    .*bool2s(abs(zv-zc)>r0);
draw_flow_conduction(x,y,Phi,([-1.5:0.25:1.5]*abs(r0*Ex0)-conj(Ex0)*real(zc)),Psi,[-10:0.25:10]*abs(r0*Ex0),xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
draw_points(real(zc),imag(zc));
