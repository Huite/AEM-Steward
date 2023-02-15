clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.4.1, Figure 3.21a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// type of boundary/interface condition
function [solve_type]=solve_type_Phi_uniform()
  solve_type="solution type: potential uniform"; 
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

function slit_flow_specify(z1,z2,solve_type,solve_value,N,M)
// z1,z2: endpoints of the slit element
// solve_type: type of boundary condition
// solve_value: boundary condition takes on values associated with
//   solve_type==solve_type_Phi_uniform()    set to [] since unused
// N: number of strength coefficients
// M: number of control points where boundary conditions are applied
  global slit;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(3*N/2); end
  sliti.z1=z1;
  sliti.z2=z2;
  sliti.L=abs(z2-z1)/2;
  sliti.expminusialpha=conj(z2-z1)/abs(z2-z1);
  sliti.solve_type=solve_type;
  sliti.solve_theta=%pi*([1:M]'-0.5)/M;
  sliti.solve_Zp=exp( %i*sliti.solve_theta);
  sliti.solve_Zm=exp(-%i*sliti.solve_theta);
  sliti.solve_zeta=0.5*(sliti.solve_Zp+sliti.solve_Zp.^(-1));
  sliti.solve_z=(sliti.z2+sliti.z1)/2+sliti.solve_zeta*(sliti.z2-sliti.z1)/2;
  // setup matrices to be used in solving system of equations to determine coefficients
  if solve_type==solve_type_Phi_uniform() then
    sliti.A=ones(sliti.solve_theta);
    for n=1:N sliti.A=[sliti.A,cos(n*sliti.solve_theta)]; end
  end
  sliti.coef=zeros(N,1);
  sliti.coef_Q=0; 
  sliti.coef_Gamma=0; 
  if max(size(slit))==0 then
    slit=list(sliti);
  else
    slit($+1)=sliti;
  end
endfunction

function [omega]=slit_flow_omega(z,excluded)
  global slit;
  if ~exists('excluded') then excluded=[]; end
  omega=zeros(z);
  for ielt=1:max(size(slit))
    if find(ielt==excluded)==[] then 
      zeta=(z-(slit(ielt).z2+slit(ielt).z1)/2)*2/(slit(ielt).z2-slit(ielt).z1);
      Z=zeta+((zeta+1).^(0.5)).*((zeta-1).^(0.5));
      logZ=log(Z);
      omega=omega+slit(ielt).coef_Q*(logZ/(2*%pi)-1);
      omega=omega+slit(ielt).coef_Gamma*%i*(logZ/(2*%pi)-1);
      for n=1:max(size(slit(ielt).coef))
        omega=omega+slit(ielt).coef(n)*Z.^(-n);
      end
    end
  end
endfunction

function [v]=slit_flow_v(z,excluded)
  global slit;
  if ~exists('excluded') then excluded=[]; end
  v=zeros(z);
  for ielt=1:max(size(slit))
    if find(ielt==excluded)==[] then
      zeta=(z-(slit(ielt).z2+slit(ielt).z1)/2)*2/(slit(ielt).z2-slit(ielt).z1);
      zetaplus1 =zeta+1;  zetaplus1 =zetaplus1 +%i*1e-6*(abs(zetaplus1) <1e-6);
      zetaminus1=zeta-1;  zetaminus1=zetaminus1+%i*1e-6*(abs(zetaminus1)<1e-6);
      Z=zeta+((zeta+1).^(0.5)).*((zeta-1).^(0.5));
      v=v-slit(ielt).coef_Q*conj(1/%pi*((zetaplus1).^(-.5)).*((zetaminus1).^(-.5))./(slit(ielt).z2-slit(ielt).z1));
      v=v+slit(ielt).coef_Gamma*%i*conj(1/%pi*((zetaplus1).^(-.5)).*((zetaminus1).^(-.5))./(slit(ielt).z2-slit(ielt).z1));
      for n=1:max(size(slit(ielt).coef))
        v=v+conj(slit(ielt).coef(n)*2*n/(slit(ielt).z2-slit(ielt).z1)*Z.^(-n).*((zetaplus1).^(-.5)).*((zetaminus1).^(-.5)));
      end
    end
  end
endfunction

function [omega]=omega_add_slit_flow(z,excluded)
  omega=uniform_flow_omega(z)...
       +slit_flow_omega(z,excluded);
endfunction

function [change]=slit_flow_solve()
  global slit
  change=0;
  for ielt=1:max(size(slit))
    omega_beforesolve=slit_flow_omega(slit(ielt).solve_z,[1:ielt-1,ielt+1:max(size(slit))]);
    M=max(size(slit(ielt).solve_z));
    if slit(ielt).solve_type==solve_type_Phi_uniform() then
      Phi_solve_z=real(omega_add_slit_flow(slit(ielt).solve_z,[ielt]));
      x=2/M*slit(ielt).A'*(-Phi_solve_z);
      slit(ielt).coef_Phi=-x(1)/2;   
      slit(ielt).coef=x(2:$);
    end
    omega_aftersolve=slit_flow_omega(slit(ielt).solve_z,[1:ielt-1,ielt+1:max(size(slit))]);
    change=max(change,max(abs(real(omega_aftersolve-omega_beforesolve))));
  end
endfunction

function slit_flow_error()
  global slit
  if max(size(slit))<1 then return end
  for ielt=1:size(slit)
    z=slit(ielt).solve_z;
    omega=uniform_flow_omega(z)+slit_flow_omega(z);
    v    =uniform_flow_v(z)    +slit_flow_v(z);
    if slit(ielt).solve_type==solve_type_Phi_uniform() then
      slit(ielt).residual_error=real(omega)-slit(ielt).coef_Phi;
    end
  end
  net_error=[];
  for ielt=1:size(slit)
    error_rmse_ielt = sqrt( sum(    ( slit(ielt).residual_error ).^2 )/max(size( slit(ielt).residual_error )) );     
    error_abs_ielt  =       sum( abs( slit(ielt).residual_error )    )/max(size( slit(ielt).residual_error ))  ;
    net_error=[net_error;slit(ielt).residual_error];
    printf("Slit %3d error rmse: %10.4e abs: %10.4e %s\n",ielt,error_rmse_ielt,error_abs_ielt,slit(ielt).solve_type);
  end
  if max(size(net_error))>0 then printf("Slit net error rmse: %10.4e average absolute value: %10.4e \n",...
     sqrt( sum(net_error.^2)/max(size(net_error)) ), sum(abs(net_error))/max(size(net_error)) ); end
endfunction

function slit_flow_draw(draw_geometry)
  global slit;
  if ~exists('draw_geometry') then draw_geometry=%F; end
  for ielt=1:max(size(slit))
    draw_segments(real(slit(ielt).z1),imag(slit(ielt).z1),real(slit(ielt).z2),imag(slit(ielt).z2));
    if draw_geometry then
      draw_points(real([slit(ielt).z1;slit(ielt).z2]),imag([slit(ielt).z1;slit(ielt).z2]));
      draw_points(real(slit(ielt).solve_z),imag(slit(ielt).solve_z));
    end
  end
endfunction

xmin=0; xmax=100; ymin=0; ymax=100;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

// set the uniform flow components for Phi0 and v0 instead of solving to match boundary conditions
uniform_flow_specify([],[]); global uniform; uniform.coef_Phi0=.9; uniform.coef_v0=.008;

solve_type=solve_type_Phi_uniform();   solve_value=[];
z1=30+%i*15; z2=70+%i*85;
slit_flow_specify(z1,z2,solve_type,solve_value);

printf("Solving for unknown coefficients\n");
count=0; change=100;
while count<5000 & change>1e-12 do
  count=count+1;
  change=0;
  change=max(change,slit_flow_solve());
  printf("Solving iteration %3d change %12.5e\n",count,change);
end
slit_flow_error();

omega=uniform_flow_omega(z)+slit_flow_omega(z);
v    =uniform_flow_v(zv)   +slit_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[0:.025:1],imag(omega),[-10:.025:10],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
slit_flow_draw();
