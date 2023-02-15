clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.3.2, Figure 1.51a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=-20; xmax=0; ymin=-10; ymax=10;
xrange=xmin+(xmax-xmin)*[0:1/400:1]; yrange=ymin+(ymax-ymin)*[0:1/400:1]'; x=(xrange'*ones(yrange'))'; y=(ones(xrange' )*yrange')';
phi0=0.1; theta0=45*%pi/180; k=2*%pi/5; R=1.0;
phi=phi0*(exp(%i*k*(x*cos(theta0)+y*sin(theta0)))+R*exp(%i*k*(-x*cos(theta0)+y*sin(theta0))));

draw_waves(xrange,yrange,phi,[],%F);
xtitle(title_author,"x","y");
