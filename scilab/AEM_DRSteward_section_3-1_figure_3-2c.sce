clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.1, Figure 3.2c";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=-1; xmax=1; ymin=-1; ymax=1;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

v0=1*exp(%i*0*%pi);
omega=-conj(v0)*z;
v=v0*ones(zv);
Gamma=1; zp=.6-%i*.6;
omega=omega+%i*Gamma/(2*%pi)*log(z-zp);
v    =v+%i*Gamma/(2*%pi)*conj((zv-zp).^(-1));
Gamma=1; zp=.4+%i*.6;
omega=omega+%i*Gamma/(2*%pi)*log(z-zp);
v    =v+%i*Gamma/(2*%pi)*conj((zv-zp).^(-1));
Gamma=1; zp=-.4-%i*.4;
omega=omega+%i*Gamma/(2*%pi)*log(z-zp);
v    =v+%i*Gamma/(2*%pi)*conj((zv-zp).^(-1));

draw_flow_conduction(x,y,real(omega),[-1.5:.1:2.5],imag(omega),[-10:.1:10],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
