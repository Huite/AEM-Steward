clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 2.3.2, Figure 2.14b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=0; xmax=1; ymin=0; ymax=1;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/30:1/15:1]; yv=ymin+(ymax-ymin)*[1/30:1/15:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';
omega=z.^2;
v    =-2*conj(zv.^(2-1));
draw_flow_conduction(x,y,real(omega),[-1:.1:1],imag(omega),[-5:.1:5],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
