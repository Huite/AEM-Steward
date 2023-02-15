clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.2.2, Figure 3.11a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// type of boundary/interface condition
function [solve_type]=solve_type_Phi_specified()
  solve_type="solution type: potential specified"; 
endfunction

function circle_flow_specify(zc,r0,solve_type,solve_value,N,M)
// zc: center of circle
// r0: radius of circle
// solve_type: type of boundary condition
// solve_value: boundary condition takes on values associated with
//   solve_type==solve_type_Phi_specified():    potential
// N: number of strength coefficients
// M: number of control points where boundary conditions are applied
  global circle;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(3*N); end
  circlei.zc=zc;
  circlei.r0=r0;
  circlei.solve_type=solve_type;
  circlei.N =N;
  circlei.M =M;
  circlei.solve_theta=2*%pi*[0:1/M:1-1/M]';
  circlei.solve_Z=exp(%i*circlei.solve_theta);
  circlei.solve_z=circlei.zc+circlei.r0*exp(%i*circlei.solve_theta);
  if solve_type==solve_type_Phi_specified() then
    circlei.solve_Phi=solve_value;
    circlei.coef_in =zeros(N+1,1);
  end
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
      omegam=zeros(z);
      for n=0:max(size(circle(ielt).coef_in))-1
        omegam=omegam+circle(ielt).coef_in(n+1)*Z.^(n);
      end
      omega=omega+omegam.*(r<1);
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
      vm=zeros(z);
      for n=1:max(size(circle(ielt).coef_in))-1
        vm=vm-conj((n)/circle(ielt).r0*circle(ielt).coef_in(n+1)*Z.^(n-1));
      end
      v=v+vm.*(r<1);
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

xmin=-4; xmax=4; ymin=-4; ymax=4;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/32:1/16:1]; yv=ymin+(ymax-ymin)*[1/32:1/16:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

zc=0; r0=4; Phi_circle=0;
circle_flow_specify(zc,r0,solve_type_Phi_specified(),Phi_circle);  

// functions for other point-sinks inside circle
zp=2*sqrt(2)*exp(%i*2*%pi*[1/8:1/4:1]);  Q=-1;
function omega=other_omega(z);
  omega=zeros(z);
  for ip=1:max(size(zp))
    omega=omega+Q/(2*%pi)*log(z-zp(ip));
  end
endfunction
function v=other_v(zv);
  v=zeros(zv);
  for ip=1:max(size(zp))
    v=v-conj(Q/(2*%pi)*(zv-zp(ip)).^(-1));
  end
endfunction

// solve for coefficients inside the circle
global circle
jcirc=1;  M=max(size(circle(jcirc).solve_z));
b=circle(jcirc).solve_Phi-real(other_omega(circle(jcirc).solve_z));
circle(jcirc).coef_in(1)=sum(b)/M;
for n=1:max(size(circle(jcirc).coef_in))-1;
  circle(jcirc).coef_in(n+1)=2/M*exp(%i*n*circle(jcirc).solve_theta')*b;
end
// residual error
tolerance=1e-8; ielt=1;
z_minus=circle(ielt).zc+(1-tolerance)*circle(ielt).r0*exp(%i*circle(ielt).solve_theta);
omega_minus=circle_flow_omega(z_minus)+other_omega(z_minus);
circle(ielt).residual_error=real(omega_minus)-circle(ielt).solve_Phi;
error_rmse_ielt = sqrt( sum(    ( circle(ielt).residual_error ).^2 )/max(size( circle(ielt).residual_error )) );     
error_abs_ielt  =       sum( abs( circle(ielt).residual_error )    )/max(size( circle(ielt).residual_error ))  ;
printf("Circle %3d error rmse: %10.4e abs: %10.4e %s\n",ielt,error_rmse_ielt,error_abs_ielt,circle(ielt).solve_type);

omega=circle_flow_omega(z)+other_omega(z).*(abs(z-zc)<=r0);
v    =circle_flow_v(zv)   +other_v(zv).*(abs(zv-zc)<=r0);

draw_flow_conduction(x,y,real(omega),[0:.05:.6],imag(omega),[-10:.05:10],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
circle_flow_draw(%F);
draw_points(real(zp),imag(zp));
