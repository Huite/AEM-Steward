clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.2.1, Figure 3.8";
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
  circlei.coef_out=zeros(N  ,1);
  circlei.coef_Q  =0;
  if max(size(circle))==0 then
    circle=list(circlei);
  else
    circle($+1)=circlei;
  end
endfunction

function [omega]=circle_flow_omega(z,excluded)
  global circle;
  if ~exists('excluded') then excluded=[]; end
  omega=zeros(z);
  for ielt=1:max(size(circle))
    if find(ielt==excluded)==[] then 
      Z=(z-circle(ielt).zc)/circle(ielt).r0; r=abs(Z); Z=Z+1e-6*(r<1e-6);
      omegap=circle(ielt).coef_Q*(log(Z)/(2*%pi)-1);
      for n=1:max(size(circle(ielt).coef_out))
        omegap=omegap+circle(ielt).coef_out(n)*Z.^(-n);
      end
      omega=omega+omegap.*(r>=1);
    end
  end
endfunction

function [v]=circle_flow_v(z,excluded)
  global circle;
  if ~exists('excluded') then excluded=[]; end
  v=zeros(z);
  for ielt=1:max(size(circle))
    if find(ielt==excluded)==[] then
      Z=(z-circle(ielt).zc)/circle(ielt).r0; r=abs(Z); Z=Z+1e-6*(r<1e-6);
      vp=-circle(ielt).coef_Q*1/(2*%pi*circle(ielt).r0)*conj(Z.^(-1));
      for n=1:max(size(circle(ielt).coef_out))
        vp=vp-conj(-n/circle(ielt).r0*circle(ielt).coef_out(n)*Z.^(-n-1));
      end
      v=v+vp.*(r>=1);
    end
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

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/36:1/18:1]; yv=ymin+(ymax-ymin)*[1/36:1/18:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

zc=0; r0=1.0; Phi_circle=0;
circle_flow_specify(zc,r0);  

// uncomment to display influence functions for series outside element
global circle
//circle(1).coef_Q=1;
circle(1).coef_out(1)=1;
//circle(1).coef_out(1)=%i;
//circle(1).coef_out(2)=1;
//circle(1).coef_out(2)=%i;
//circle(1).coef_out(3)=1;
//circle(1).coef_out(3)=%i;

omega=circle_flow_omega(z);
v    =circle_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[-1:.1:1],imag(omega),[-1.5:.1:1.5],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
circle_flow_draw();
