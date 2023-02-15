clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.4, Figure 3.18b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/24:1/12:1]; yv=ymin+(ymax-ymin)*[1/24:1/12:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

v0=1*exp(%i*%pi/3);
omega=-0.5*conj(v0)*( (((z+1).^(0.5)).*((z-1).^(0.5))) + z)...
      -0.5*v0      *( (((z+1).^(0.5)).*((z-1).^(0.5))) + z).^(-1) ;
v    = 0.5*v0      *conj(  (((zv+1).^(0.5)).*((zv-1).^(0.5)) + zv)./(((zv+1).^(0.5)).*((zv-1).^(0.5)))  )...
      -0.5*conj(v0)*conj( ((((zv+1).^(0.5)).*((zv-1).^(0.5)) + zv).*(((zv+1).^(0.5)).*((zv-1).^(0.5)))).^(-1) );

draw_flow_conduction(x,y,real(omega),[-2.2:.2:2]*abs(v0),imag(omega),[-10:0.2:10]*abs(v0),xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
draw_segments(-1,0,1,0);
