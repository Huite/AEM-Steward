clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.4, Figure 3.19a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/24:1/12:1]; yv=ymin+(ymax-ymin)*[1/24:1/12:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

Q=1;
omega=Q/(2*%pi)*log( (((z+1).^(0.5)) .* ((z-1).^(0.5))) + z);
v    =-Q/(2*%pi)*( conj(  (zv+1).^(-0.5) .* (zv-1).^(-0.5) ));
v=v.*(~(abs(imag(zv))<0.2&real(zv)<1.2&real(zv)>-1.2));

draw_flow_conduction(x,y,real(omega),[0:.025:0.25]*Q,imag(omega),[-5:0.025:5]*Q,xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
draw_segments(-1,0,1,0);
