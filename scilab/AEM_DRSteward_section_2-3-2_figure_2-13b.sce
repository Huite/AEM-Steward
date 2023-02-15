clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 2.3.2, Figure 2.13b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

z1=2*exp(%i*%pi/12);
z2=5*exp(%i*%pi/6);
xmin=0; xmax=10; ymin=0; ymax=10;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/30:1/15:1]; yv=ymin+(ymax-ymin)*[1/30:1/15:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';
draw_flow_conduction(x,y,[]);
xtitle(title_author,"x","y");

printf("z_1     = %.4f e^{ i %.4f}\n",abs(z1),atan(imag(z1),real(z1)));
printf("z_2     = %.4f e^{ i %.4f}\n",abs(z2),atan(imag(z2),real(z2)));
printf("z_1 z_2 = %.4f e^{ i %.4f}\n",abs(z1*z2),atan(imag(z1*z2),real(z1*z2)));
draw_points(real([z1;z2;z1*z2]),imag([z1;z2;z1*z2]));
