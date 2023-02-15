clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.2.3, Figure 1.26b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/24:1/12:1]; yv=ymin+(ymax-ymin)*[1/24:1/12:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';
z0=0; r0=1; v0=1; vx0=-v0; Gamma=6;
omega=-v0*z*exp(-%i*%pi)-v0*exp(%i*%pi)*r0^2*(z-z0).^(-1);
omega=omega+%i*Gamma/(2*%pi)*log(z);
vx=vx0-vx0*(real(zv).^2-imag(zv).^2)./((real(zv).^2+imag(zv).^2).^2)-Gamma/(2*%pi)*imag(zv)./(real(zv).^2+imag(zv).^2);
vy=   -vx0*(2*real(zv).*imag(zv))./((real(zv).^2+imag(zv).^2).^2)   +Gamma/(2*%pi)*real(zv)./(real(zv).^2+imag(zv).^2);
omega=omega.*bool2s(abs(z -z0)>r0)+100*bool2s(abs(z -z0)<r0);
vx   =vx   .*bool2s(abs(zv-z0)>r0);
vy   =vy   .*bool2s(abs(zv-z0)>r0);
draw_flow_conduction(x,y,real(omega),[-6:.25:4],imag(omega),[-10:.25:10],xv,yv,vx,vy);
xtitle(title_author,"x","y");
