clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.2, Figure 3.7a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

function circle_flow_specify(zc,r0,N,M)
// zc: center of circle
// r0: radius of circle
// N: number of strength coefficients
// M: number of control points where boundary conditions are applied
  global circle;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(3*N); end
  circlei.zc=zc;
  circlei.r0=r0;
  circlei.N =N;
  circlei.M =M;
  circlei.solve_theta=2*%pi*[0:1/M:1-1/M]';
  circlei.solve_Z=exp(%i*circlei.solve_theta);
  circlei.solve_z=circlei.zc+circlei.r0*exp(%i*circlei.solve_theta);
  if max(size(circle))==0 then
    circle=list(circlei);
  else
    circle($+1)=circlei;
  end
endfunction

function circle_flow_draw(draw_geometry)
  global circle;
  if ~exists('draw_geometry') then draw_geometry=%F; end
  for ielt=1:max(size(circle))
    boundary=circle(ielt).zc+circle(ielt).r0*exp(%i*[0:%pi/36:2*%pi]);
    draw_segments(real(boundary(1:$-1)),imag(boundary(1:$-1)),real(boundary(2:$)),imag(boundary(2:$)));
    if draw_geometry then
      draw_points(real(circle(ielt).zc),imag(circle(ielt).zc));
      draw_points(real(circle(ielt).solve_z),imag(circle(ielt).solve_z));
    end
  end
endfunction

xmin=0; xmax=100; ymin=0; ymax=100;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/32:1/16:1]; yv=ymin+(ymax-ymin)*[1/32:1/16:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

zc=50+%i*50; r0=40;
circle_flow_specify(zc,r0);  

draw_flow_conduction(x,y,[],[],[],[],xv,yv,[],[]);
xtitle(title_author,"x","y");
circle_flow_draw(%T);
