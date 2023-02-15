clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.3, Figure 3.14a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

function ellipse_flow_specify(zc,Lmajor,Lminor,angle,N,M)
// zc: center of ellipse
// Lmajor,Lminor: 1/2 length of major and minor axes of ellipse
// angle: orientation of ellipse with respect to x-axis in units of radians
// N: number of strength coefficients
// M: number of control points where boundary conditions are applied
  global ellipse;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(3*N); end
  ellipsei.zc=zc;
  ellipsei.Lmajor=Lmajor;
  ellipsei.Lminor=Lminor;
  ellipsei.Lfoci=sqrt(ellipsei.Lmajor^2-ellipsei.Lminor^2);
  ellipsei.F=ellipsei.Lfoci/(ellipsei.Lmajor+ellipsei.Lminor);
  ellipsei.alpha=angle;
  ellipsei.solve_theta=2*%pi*[0:1/M:1-1/M]';
  ellipsei.solve_Z=exp(%i*2*%pi*[0:1/M:1-1/M]');
  ellipsei.solve_zeta=0.5*(ellipsei.solve_Z+ellipsei.F^2*ellipsei.solve_Z.^(-1));
  ellipsei.solve_z=ellipsei.zc+ellipsei.solve_zeta*(ellipsei.Lmajor+ellipsei.Lminor)*exp(%i*ellipsei.alpha);
  if max(size(ellipse))==0 then
    ellipse=list(ellipsei);
  else
    ellipse($+1)=ellipsei;
  end
endfunction

function ellipse_flow_draw(draw_geometry)
  global ellipse;
  if ~exists('draw_geometry') then draw_geometry=%F; end
  for ielt=1:max(size(ellipse))
    drawZ=exp(%i*[0:%pi/36:2*%pi]);
    drawzeta=0.5*(drawZ+ellipse(ielt).F^2*drawZ.^(-1));
    drawz=ellipse(ielt).zc+drawzeta*(ellipse(ielt).Lmajor+ellipse(ielt).Lminor)*exp(%i*ellipse(ielt).alpha);
    draw_segments(real(drawz(1:$-1)),imag(drawz(1:$-1)),real(drawz(2:$)),imag(drawz(2:$)));
    if draw_geometry then
      zdm=ellipse(ielt).zc-ellipse(ielt).Lfoci*exp(%i*ellipse(ielt).alpha);
      zdp=ellipse(ielt).zc+ellipse(ielt).Lfoci*exp(%i*ellipse(ielt).alpha);
      draw_points(real([ellipse(ielt).zc,zdm,zdp]),imag([ellipse(ielt).zc,zdm,zdp]));
      draw_points(real(ellipse(ielt).solve_z),imag(ellipse(ielt).solve_z));
     end
  end
endfunction

xmin=0; xmax=100; ymin=0; ymax=100;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

zc=50+%i*50; Lmajor=40; Lminor=24; angle=%pi/3;
ellipse_flow_specify(zc,Lmajor,Lminor,angle);

draw_flow_conduction(x,y,[],[],[],[],xv,yv,[],[]);
xtitle(title_author,"x","y");
ellipse_flow_draw(%T);
