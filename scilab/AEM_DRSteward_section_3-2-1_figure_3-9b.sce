clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.2.1, Figure 3.9b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// type of boundary/interface condition
function [solve_type]=solve_type_Phi_specified()
  solve_type="solution type: potential specified"; 
endfunction

function uniform_flow_specify(solve_z,solve_Phi)
// solve_z: boundary condtion column vector with location of control points
// solve_Phi: boundary condition column vector with specified potential at the control points
// uniform.coef_Phi0: coefficient Phi0 in uniform potential function
// uniform.coef_v0: coefficient v0 in uniform vector field
  global uniform
  uniform.solve_z=solve_z;
  uniform.solve_Phi=solve_Phi;
  uniform.coef_Phi0=0;
  uniform.coef_v0=0;
endfunction
  
function [omega]=uniform_flow_omega(z)
  global uniform
  omega=-conj(uniform.coef_v0)*z+uniform.coef_Phi0;
endfunction

function [v]=uniform_flow_v(z)
  global uniform
  v=uniform.coef_v0*ones(z);
endfunction

function [omega]=omega_add_uniform_flow(z)
  omega=circle_flow_omega(z);
endfunction

function [change]=uniform_flow_solve()
  global uniform
  omega_beforesolve=uniform_flow_omega(uniform.solve_z);
  A=[ones(uniform.solve_z),-real(uniform.solve_z),-imag(uniform.solve_z)]; 
  b=real(uniform.solve_Phi-omega_add_uniform_flow(uniform.solve_z));
  coef=A\b;
  uniform.coef_Phi0=coef(1);
  uniform.coef_v0=coef(2)+%i*coef(3);
  omega_aftersolve=uniform_flow_omega(uniform.solve_z);
  change=max(abs(real(omega_aftersolve-omega_beforesolve)));
endfunction

function uniform_flow_error()
  global uniform
  z=uniform.solve_z;
  omega=uniform_flow_omega(z)+circle_flow_omega(z);
  residual=uniform.solve_Phi-real(omega);
  printf("Uniform flow error: rmse: %11.5e average abs: %11.5e\n",...
         sqrt( sum(     residual.^2 )/max(size(residual)) ) ,...
               sum( abs(residual)   )/max(size(residual))   );
endfunction

function uniform_flow_draw()
  global uniform;
  if max(size(uniform.solve_z))>0 then
    draw_points(real(uniform.solve_z),imag(uniform.solve_z));
  end
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
    circlei.coef_out=zeros(N  ,1);
    circlei.coef_Q  =0;
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

function [omega]=omega_add_circle_flow(z,excluded)
  omega=uniform_flow_omega(z)...
       +circle_flow_omega(z,excluded);
endfunction

function [change]=circle_flow_solve()
  global circle
  change=0;
  for ielt=1:max(size(circle))
    omega_beforesolve=circle_flow_omega(circle(ielt).solve_z,[1:ielt-1,ielt+1:max(size(circle))]);
    if circle(ielt).solve_type==solve_type_Phi_specified() then
      Phi_circle=circle(ielt).solve_Phi;
      Phi_solve_z=real(omega_add_circle_flow(circle(ielt).solve_z,[ielt]));
      M=max(size(circle(ielt).solve_z));
      circle(ielt).coef_Q=-sum(Phi_circle-Phi_solve_z)/M;
      circle(ielt).coef_out=zeros(circle(ielt).coef_out);
      for m=1:M
        for n=1:max(size(circle(ielt).coef_out))
          circle(ielt).coef_out(n)=circle(ielt).coef_out(n)+2/M*(Phi_circle-Phi_solve_z(m))*exp(%i*(n)*circle(ielt).solve_theta(m));
        end
      end
    end
    omega_aftersolve=circle_flow_omega(circle(ielt).solve_z,[1:ielt-1,ielt+1:max(size(circle))]);
    change=max(change,max(abs(real(omega_aftersolve-omega_beforesolve))));
  end
endfunction

function circle_flow_error()
  global circle
  if max(size(circle))<1 then return end
  for ielt=1:size(circle)
	tolerance=1e-8;
    z_plus=circle(ielt).zc+(1+tolerance)*circle(ielt).r0*exp(%i*circle(ielt).solve_theta);
    omega_plus=uniform_flow_omega(z_plus)+circle_flow_omega(z_plus);
    if circle(ielt).solve_type==solve_type_Phi_specified() then
      circle(ielt).residual_error=real(omega_plus)-circle(ielt).solve_Phi;
    end
  end
  net_error=[];
  for ielt=1:size(circle)
    error_rmse_ielt = sqrt( sum(    ( circle(ielt).residual_error ).^2 )/max(size( circle(ielt).residual_error )) );     
    error_abs_ielt  =       sum( abs( circle(ielt).residual_error )    )/max(size( circle(ielt).residual_error ))  ;
    net_error=[net_error;circle(ielt).residual_error];
    printf("Circle %3d error rmse: %10.4e abs: %10.4e %s\n",ielt,error_rmse_ielt,error_abs_ielt,circle(ielt).solve_type);
  end
  if max(size(net_error))>0 then printf("Circle net error rmse: %10.4e average absolute value: %10.4e \n",...
     sqrt( sum(net_error.^2)/max(size(net_error)) ), sum(abs(net_error))/max(size(net_error)) ); end
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

xmin=-1; xmax=1; ymin=-1; ymax=1;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/32:1/16:1]; yv=ymin+(ymax-ymin)*[1/32:1/16:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

solve_z=[0]; 
solve_Phi=[0]';
uniform_flow_specify(solve_z,solve_Phi);
zc=-0.5; r0=0.2; Phi_circle=-0.5;
circle_flow_specify(zc,r0,solve_type_Phi_specified(),Phi_circle);  
zc=0.5; r0=0.2; Phi_circle=0.5;
circle_flow_specify(zc,r0,solve_type_Phi_specified(),Phi_circle);  

printf("Solving for unknown coefficients\n");
count=0; change=100;
while count<5000 & change>1e-12 do
  count=count+1;
  change=uniform_flow_solve();
  change=max(change,circle_flow_solve());
  printf("Solving iteration %3d change %12.5e\n",count,change);
end
uniform_flow_error();
circle_flow_error();

omega=uniform_flow_omega(z)+circle_flow_omega(z);
v    =uniform_flow_v(zv)   +circle_flow_v(zv);
global circle
for icirc=1:size(circle)
  omega=omega+(circle(icirc).solve_Phi-omega).*(abs(z-circle(icirc).zc)<circle(icirc).r0);
  v=v+(-v).*(abs(zv-circle(icirc).zc)<circle(icirc).r0);
end

draw_flow_conduction(x,y,real(omega),[-0.5-1e-8:.1:0.5],imag(omega),[-10:.1:10],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
circle_flow_draw(%F);
uniform_flow_draw();
