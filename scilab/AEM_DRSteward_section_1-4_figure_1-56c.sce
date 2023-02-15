clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.4, Figure 1.56c";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=-4; xmax=4; ymin=-8; ymax=0;
xrange =xmin+(xmax-xmin)*[0:1/400:1];  yrange= ymin+(ymax-ymin)*[0:1/400:1]';  x=(xrange'*ones(yrange'))';    y=(ones(xrange' )*yrange')';
xrangev=xmin+(xmax-xmin)*[1/16:1/8:1]; yrangev=ymin+(ymax-ymin)*[1/16:1/8:1]'; xv=(xrangev'*ones(yrangev'))'; yv=(ones(xrangev')*yrangev')';

ysoil=0; ysat=-4;
density_water=1000;              // density of water kg/m^2
porosity     =0.30;              // porosity of compacted soil
density_soil =2650*(1-porosity); // density of soil kg/m^3 (2.65g/cm^3 for rocks) minus air
g            =9.81;              // acceleration of gravity m/s^2
nu           =0.3;               // Poisson ratio of sand/soil
E            =80*10^6;           // Elasticity of sand/gravel Pa=kg/(ms^2)
mu           =E/(2*(1+nu));      // Modulus of rigidity of sand/soil

// dry soil above saturated soil
sigmax    =nu/(1-nu)*density_soil*g*(y -ysoil).*(y >ysat);  sigmay    =density_soil*g*(y -ysoil).*(y >ysat);  tauxy     =zeros(x);
vec.sigmax=nu/(1-nu)*density_soil*g*(yv-ysoil).*(yv>ysat);  vec.sigmay=density_soil*g*(yv-ysoil).*(yv>ysat);  vec.tauxy =zeros(xv);
vec.ux=zeros(xv);  vec.uy=(1/(2*mu)*(1-2*nu)/(1-nu)*density_soil*g*(yv-ysoil).^2/2).*(yv>ysat);

// soil and water below saturated interface
sigmax_sat=nu/(1-nu)*density_soil*g*(ysat-ysoil);
sigmax    =sigmax    +(sigmax_sat+(nu/(1-nu)*(density_soil+(-1+porosity)*density_water)+density_water)*g*(y -ysat)).*(y <=ysat);
vec.sigmax=vec.sigmax+(sigmax_sat+(nu/(1-nu)*(density_soil+(-1+porosity)*density_water)+density_water)*g*(yv-ysat)).*(yv<=ysat);
sigmay_sat=density_soil*g*(ysat-ysoil);
sigmay    =sigmay    +(sigmay_sat+(density_soil+porosity*density_water)*g*(y -ysat)).*(y <=ysat);
vec.sigmay=vec.sigmay+(sigmay_sat+(density_soil+porosity*density_water)*g*(yv-ysat)).*(yv<=ysat);
uy_sat=(1/(2*mu)*(1-2*nu)/(1-nu)*density_soil*g*(ysat-ysoil).^2/2);
vec.uy=vec.uy+(uy_sat+1/(2*mu)*(1-2*nu)/(1-nu)*(density_soil+(-1+porosity)*density_water)*g*(yv-ysat).^2/2).*(yv<=ysat);
  
// visualize stress
draw_stress(xrange,yrange,sigmax,sigmay,tauxy,[-200000:10000:0],xrangev,yrangev,vec,%T);
// visualize displacement
//draw_stress(xrange,yrange,sigmax,sigmay,tauxy,[],xrangev,yrangev,vec,%F);

draw_segments([xmin;xmin],[ysoil;ysat],[xmax;xmax],[ysoil;ysat]);
xtitle(title_author,"x","y");
