clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.4.3, Figure 3.27b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// types of boundary/interface condition
function [solve_type]=solve_type_Psi_specified()
  solve_type="solution type: stream function specified"; 
endfunction
function [solve_type]=solve_type_Kutta()
  solve_type="solution type: Kutta condition and stream function uniform"; 
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
//   solve_type==solve_type_Psi_specified()  stream function
//   solve_type==solve_type_Kutta()          %T: first endpoint, %F: second endpoint (default %T if [])
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
  if sliti.solve_type==solve_type_Psi_specified() then
    sliti.solve_Psi=solve_value;
  elseif sliti.solve_type==solve_type_Kutta() then
    if solve_value==[] then solve_value=%T; end
    sliti.solve_Kutta_end1=solve_value;
  end
  sliti.solve_theta=%pi*([1:M]'-0.5)/M;
  sliti.solve_Zp=exp( %i*sliti.solve_theta);
  sliti.solve_Zm=exp(-%i*sliti.solve_theta);
  sliti.solve_zeta=0.5*(sliti.solve_Zp+sliti.solve_Zp.^(-1));
  sliti.solve_z=(sliti.z2+sliti.z1)/2+sliti.solve_zeta*(sliti.z2-sliti.z1)/2;
  // setup matrices to be used in solving system of equations to determine coefficients
  if  solve_type==solve_type_Psi_specified() then
    sliti.A=ones(sliti.solve_theta);
    for n=1:N sliti.A=[sliti.A,cos(n*sliti.solve_theta)]; end
  elseif solve_type==solve_type_Kutta() then
    sliti.A=[];
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

function [v]=v_add_slit_flow(z,excluded)
  v=uniform_flow_v(z)...
   +slit_flow_v(z,excluded);
endfunction

function [change]=slit_flow_solve()
  global slit
  change=0;
  for ielt=1:max(size(slit))
    omega_beforesolve=slit_flow_omega(slit(ielt).solve_z,[1:ielt-1,ielt+1:max(size(slit))]);
    M=max(size(slit(ielt).solve_z));
    if  slit(ielt).solve_type==solve_type_Psi_specified() then
      Psi_solve_z=imag(omega_add_slit_flow(slit(ielt).solve_z,[ielt]));
      x=2/M*slit(ielt).A'*(slit(ielt).solve_Psi-Psi_solve_z);
      slit(ielt).coef_Gamma=-x(1)/2;
      slit(ielt).coef=%i*x(2:$);      
    elseif slit(ielt).solve_type==solve_type_Kutta() then
      Psi_solve_z=imag(omega_add_slit_flow(slit(ielt).solve_z,[ielt]));
      x=2/M*slit(ielt).A'*(-Psi_solve_z);
      slit(ielt).coef=%i*x;
      N=max(size(slit(ielt).coef));
      if slit(ielt).solve_Kutta_end1 then
        slit(ielt).coef_Gamma=2*%pi*sum([1:N]'.*((-1).^[1:N]').*imag(slit(ielt).coef));
      else
        slit(ielt).coef_Gamma=2*%pi*sum([1:N]'.*imag(slit(ielt).coef));
      end
      slit(ielt).Psiav=sum(Psi_solve_z)/max(size(Psi_solve_z))-slit(ielt).coef_Gamma;
      // code to adjust the specified stream function of other slits to be that of this element
      adjust_Psi=%T;
      if adjust_Psi then
        solve_Psi=imag(slit(ielt).A*slit(ielt).coef)+Psi_solve_z-slit(ielt).coef_Gamma;
        solve_Psiav=sum(solve_Psi)/max(size(solve_Psi));
		jslitelt=[1:ielt-1,ielt+1:size(slit)];
        for jslit=jslitelt
          slit(jslit).solve_Psi=solve_Psiav;
        end
        printf("adjusting solve_Psi: %f ",solve_Psiav);
      end
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
    if slit(ielt).solve_type==solve_type_Kutta() then
      slit(ielt).residual_error=imag(omega)-slit(ielt).Psiav;
    elseif slit(ielt).solve_type==solve_type_Psi_specified() then
      slit(ielt).residual_error=imag(omega)-slit(ielt).solve_Psi;
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

xmin=-1.2; xmax=1.2; ymin=-1.2; ymax=1.2;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/24:1/12:1]; yv=ymin+(ymax-ymin)*[1/24:1/12:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

// set the uniform flow components for Phi0 and v0 instead of solving to match boundary conditions
uniform_flow_specify([],[]); global uniform; uniform.coef_Phi0=0; uniform.coef_v0=1*exp(%i*%pi);

// locate slit elements along the geometry of the wing
L=0.75; alpha=9.5/180*%pi; N=20;
zwing=L*(-%i/tan(alpha)+1/sin(alpha)*exp(%i*(%pi/2+[alpha:-alpha*2/N:-alpha]')))*exp(%i*alpha); 

  //  L=0.75; alpha=9.5/180*%pi; zwing=[-L;L]*exp(%i*alpha/180*%pi);

// specify left end of the element to apply the Kutta condition
solve_type=solve_type_Kutta(); solve_value=%T;
ielt=1;  slit_flow_specify(zwing(ielt),zwing(ielt+1),solve_type,solve_value);
// specify other elements to use a boundary condition with the same value of streamfunction as the leftmost element
for ielt=2:max(size(zwing))-1
  solve_type=solve_type_Psi_specified(); solve_value=3;
  slit_flow_specify(zwing(ielt),zwing(ielt+1),solve_type,solve_value);
end

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

global slit
if size(slit)>1 then Psi_wing=slit(2).solve_Psi; else Psi_wing=0; end

draw_flow_conduction(x,y,[],[],imag(omega),[-2+Psi_wing:0.05:2],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
slit_flow_draw();
