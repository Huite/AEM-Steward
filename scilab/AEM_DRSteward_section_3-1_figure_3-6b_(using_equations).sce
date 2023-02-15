clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.1, Figure 3.6b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=0; xmax=100; ymin=0; ymax=100;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

global uniform
uniform.solve_z=[10+%i*20;40+%i*80;90+%i*60];
uniform.solve_Phi=[20;14;10];
A=[ones(uniform.solve_z),-real(uniform.solve_z),-imag(uniform.solve_z)]; 
b=real(uniform.solve_Phi);
coef=A\b;
uniform.coef_Phi0=coef(1);
uniform.coef_v0=coef(2)+%i*coef(3);
omega=-conj(uniform.coef_v0)*z+uniform.coef_Phi0;
v=uniform.coef_v0*ones(zv);

draw_flow_conduction(x,y,real(omega),[7:1:22],imag(omega),[-300:1:300],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
draw_points(real(uniform.solve_z),imag(uniform.solve_z));
