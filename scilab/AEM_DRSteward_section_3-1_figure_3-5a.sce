clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.1, Figure 3.5a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

zc=0+%i*0; r0=1; v0=1*exp(%i*0*%pi); Gamma=0;
xmin=-1.5*r0+real(zc); xmax=1.5*r0+real(zc); ymin=-1.5*r0+imag(zc); ymax=1.5*r0+imag(zc);
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/24:1/12:1]; yv=ymin+(ymax-ymin)*[1/24:1/12:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

omega=-conj(v0)*z-v0*r0^2*(z-zc).^(-1)+%i*Gamma/(2*%pi)*log(z-zc);
v    =v0-conj(v0)*r0^2*(conj(zv-zc)).^(-2)+%i*Gamma/(2*%pi)*(conj(zv-zc).^(-1));
// don't draw functions inside the circle element
omega=omega.*bool2s(abs(z -zc)>=r0)+100*bool2s(abs(z-zc)<r0);
v    =v    .*bool2s(abs(zv-zc)>=r0);

draw_flow_conduction(x,y,real(omega),[-6:0.25:6],imag(omega),[-100:0.25:100]*abs(v0),xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
draw_points(real(zc),imag(zc));
