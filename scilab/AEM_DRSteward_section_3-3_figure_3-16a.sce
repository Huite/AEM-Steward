clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.3, Figure 3.16a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// type of boundary/interface condition
function [solve_type]=solve_type_Phi_cont()
  solve_type="solution type: potential continuity"; 
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

function ellipse_flow_specify(zc,Lmajor,Lminor,angle,solve_type,solve_value,N,M)
// zc: center of ellipse
// Lmajor,Lminor: 1/2 length of major and minor axes of ellipse
// angle: orientation of ellipse with respect to x-axis in units of radians
// solve_type: type of boundary condition
// solve_value: boundary condition takes on values associated with
//   solve_type==solve_type_Phi_cont():  [alpha_plus alpha_minus]
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
  ellipsei.solve_type=solve_type;
  if solve_type==solve_type_Phi_cont() then
    ellipsei.solve_alpha_plus=solve_value(1);
    ellipsei.solve_alpha_minus=solve_value(2);
  end
  ellipsei.solve_theta=2*%pi*[0:1/M:1-1/M]';
  ellipsei.solve_Z=exp(%i*2*%pi*[0:1/M:1-1/M]');
  ellipsei.solve_zeta=0.5*(ellipsei.solve_Z+ellipsei.F^2*ellipsei.solve_Z.^(-1));
  ellipsei.solve_z=ellipsei.zc+ellipsei.solve_zeta*(ellipsei.Lmajor+ellipsei.Lminor)*exp(%i*ellipsei.alpha);
  ellipsei.coef_in =zeros(N,1); 
  ellipsei.coef_out=zeros(N-1,1);
  if max(size(ellipse))==0 then
    ellipse=list(ellipsei);
  else
    ellipse($+1)=ellipsei;
  end
endfunction

function [omega]=ellipse_flow_omega(z,excluded)
  global ellipse;
  if ~exists('excluded') then excluded=[]; end
  omega=zeros(z);
  for ielt=1:max(size(ellipse))
    if find(ielt==excluded)==[] then 
      zeta=(z-ellipse(ielt).zc)/(ellipse(ielt).Lmajor+ellipse(ielt).Lminor)*exp(-%i*ellipse(ielt).alpha);
      Z=zeta+((zeta+ellipse(ielt).F).^(.5)).*((zeta-ellipse(ielt).F).^(0.5));
      r=abs(Z); Z=Z+1e-6*(r<1e-6);
      omegam=zeros(z);
      for n=0:max(size(ellipse(ielt).coef_in))-1
        omegam=omegam+ellipse(ielt).coef_in(n+1)*(Z.^(n)+ellipse(ielt).F^(2*n)*Z.^(-n));
      end
      omegap=zeros(z);
      for n=1:max(size(ellipse(ielt).coef_out))
        omegap=omegap+ellipse(ielt).coef_out(n)*Z.^(-n);
      end
      omega=omega+omegam.*(r<1)+omegap.*(r>=1);
    end
  end
endfunction

function [v]=ellipse_flow_v(z,excluded)
  global ellipse;
  if ~exists('excluded') then excluded=[]; end
  v=zeros(z);
  for ielt=1:max(size(ellipse))
    if find(ielt==excluded)==[] then
      zeta=(z-ellipse(ielt).zc)/(ellipse(ielt).Lmajor+ellipse(ielt).Lminor)*exp(-%i*ellipse(ielt).alpha);
      Z=zeta+((zeta+ellipse(ielt).F).^(.5)).*((zeta-ellipse(ielt).F).^(0.5));
      r=abs(Z); Z=Z+1e-6*(r<1e-6);
      vm=zeros(z);
      for n=1:max(size(ellipse(ielt).coef_in))-1
        vm=vm-conj(ellipse(ielt).coef_in(n+1)*n*exp(-%i*ellipse(ielt).alpha)/(ellipse(ielt).Lmajor+ellipse(ielt).Lminor)...
                   *(Z.^(n)-ellipse(ielt).F^(2*n)*Z.^(-n)).*(zeta+ellipse(ielt).F).^(-.5).*(zeta-ellipse(ielt).F).^(-.5));
      end
      vp=zeros(z);
      for n=1:max(size(ellipse(ielt).coef_out))
        vp=vp+conj(ellipse(ielt).coef_out(n)*n*exp(-%i*ellipse(ielt).alpha)/(ellipse(ielt).Lmajor+ellipse(ielt).Lminor)...
                   *Z.^(-n).*(zeta+ellipse(ielt).F).^(-.5).*(zeta-ellipse(ielt).F).^(-.5));
      end
      v=v+vm.*(r<1)+vp.*(r>=1);
    end
  end
endfunction

function [omega]=omega_add_ellipse_flow(z,excluded)
  omega=uniform_flow_omega(z)...
       +ellipse_flow_omega(z,excluded);
endfunction

function [change]=ellipse_flow_solve()
  global ellipse
  global uniform
  change=0;
  for ielt=1:max(size(ellipse))
    omega_beforesolve=ellipse_flow_omega(ellipse(ielt).solve_z,[1:ielt-1,ielt+1:max(size(ellipse))]);
    M=max(size(ellipse(ielt).solve_z));
    if ellipse(ielt).solve_type==solve_type_Phi_cont() then
      ellipse(ielt).coef_cont=zeros(ellipse(ielt).coef_in);
      phi_solve_z=real(omega_add_ellipse_flow(ellipse(ielt).solve_z,[ielt]));
      for m=1:M
        ellipse(ielt).coef_cont(1)=ellipse(ielt).coef_cont(1)...
		  - (ellipse(ielt).solve_alpha_plus-ellipse(ielt).solve_alpha_minus)/ellipse(ielt).solve_alpha_minus/M*phi_solve_z(m);
        for n=1:max(size(ellipse(ielt).coef_cont))-1
          ellipse(ielt).coef_cont(n+1)=ellipse(ielt).coef_cont(n+1)...
            + (-(ellipse(ielt).solve_alpha_plus^2-ellipse(ielt).solve_alpha_minus^2)+ellipse(ielt).F^(4*n)*(ellipse(ielt).solve_alpha_plus-ellipse(ielt).solve_alpha_minus)^2)...
             /((ellipse(ielt).solve_alpha_plus+ellipse(ielt).solve_alpha_minus)^2-ellipse(ielt).F^(4*n)*(ellipse(ielt).solve_alpha_plus-ellipse(ielt).solve_alpha_minus)^2)...
             *2/M*exp( %i*n*ellipse(ielt).solve_theta(m))*phi_solve_z(m) ...
            + ellipse(ielt).F^(2*n)*2*ellipse(ielt).solve_alpha_minus*(ellipse(ielt).solve_alpha_plus-ellipse(ielt).solve_alpha_minus)...
             /((ellipse(ielt).solve_alpha_plus+ellipse(ielt).solve_alpha_minus)^2-ellipse(ielt).F^(4*n)*(ellipse(ielt).solve_alpha_plus-ellipse(ielt).solve_alpha_minus)^2)...
             *2/M*exp(-%i*n*ellipse(ielt).solve_theta(m))*phi_solve_z(m);
        end
      end   
      ellipse(ielt).coef_in(1)=-0.5*ellipse(ielt).coef_cont(1);
      F2n=(ellipse(ielt).F).^(2*[1:(max(size(ellipse(ielt).coef_cont))-1)]');
      ellipse(ielt).coef_in(2:$)=-real(ellipse(ielt).coef_cont(2:$))./(1-F2n)+%i*imag(ellipse(ielt).coef_cont(2:$))./(1+F2n)
      ellipse(ielt).coef_out=ellipse(ielt).coef_cont(2:$);
    end
    omega_aftersolve=ellipse_flow_omega(ellipse(ielt).solve_z,[1:ielt-1,ielt+1:max(size(ellipse))]);
    change=max(change,max(abs(real(omega_aftersolve-omega_beforesolve))));
  end
endfunction

function ellipse_flow_error()
  global ellipse
  if max(size(ellipse))<1 then return end
  for ielt=1:size(ellipse)
	tolerance=1e-8;
    M=max(size(ellipse(ielt).solve_z));
    Z_plus =(1+tolerance)*exp(%i*2*%pi*[0:1/M:1-1/M]');
    Z_minus=(1-tolerance)*exp(%i*2*%pi*[0:1/M:1-1/M]');
    zeta_plus =0.5*(Z_plus +ellipse(ielt).F^2*Z_plus.^(-1));
    zeta_minus=0.5*(Z_minus+ellipse(ielt).F^2*Z_minus.^(-1));
    z_plus =ellipse(ielt).zc+zeta_plus *(ellipse(ielt).Lmajor+ellipse(ielt).Lminor)*exp(%i*ellipse(ielt).alpha);
    z_minus=ellipse(ielt).zc+zeta_minus*(ellipse(ielt).Lmajor+ellipse(ielt).Lminor)*exp(%i*ellipse(ielt).alpha);
    omega_plus =uniform_flow_omega(z_plus) +ellipse_flow_omega(z_plus);
    omega_minus=uniform_flow_omega(z_minus)+ellipse_flow_omega(z_minus);
    if ellipse(ielt).solve_type==solve_type_Phi_cont() then
      ellipse(ielt).residual_error=ellipse(ielt).solve_alpha_plus*real(omega_plus)-ellipse(ielt).solve_alpha_minus*real(omega_minus);
    end
  end
  net_error=[];
  for ielt=1:size(ellipse)
    error_rmse_ielt = sqrt( sum(    ( ellipse(ielt).residual_error ).^2 )/max(size( ellipse(ielt).residual_error )) );     
    error_abs_ielt  =       sum( abs( ellipse(ielt).residual_error )    )/max(size( ellipse(ielt).residual_error )) ;
    net_error=[net_error;ellipse(ielt).residual_error];
    printf("Ellipse %3d error rmse: %10.4e abs: %10.4e %s\n",ielt,error_rmse_ielt,error_abs_ielt,ellipse(ielt).solve_type);
  end
  if max(size(net_error))>0 then printf("Ellipse net error rmse: %10.4e average absolute value: %10.4e \n",...
     sqrt( sum(net_error.^2)/max(size(net_error)) ), sum(abs(net_error))/max(size(net_error)) ); end
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

// set the uniform flow components for Phi0 and v0 instead of solving to match boundary conditions
uniform_flow_specify([],[]); global uniform; uniform.coef_Phi0=7.5; uniform.coef_v0=.05;
zc=50+%i*50; Lmajor=40; Lminor=24; angle=%pi/3; solve_alpha_plus=1/10; solve_alpha_minus=1/2.5;
ellipse_flow_specify(zc,Lmajor,Lminor,angle,solve_type_Phi_cont(),[solve_alpha_plus solve_alpha_minus]);

printf("Solving for unknown coefficients\n");
count=0; change=100;
while count<5000 & change>1e-12 do
  count=count+1;
  change=0;
  change=max(change,ellipse_flow_solve());
  printf("Solving iteration %3d change %12.5e\n",count,change);
end
ellipse_flow_error();

omega=uniform_flow_omega(z)+ellipse_flow_omega(z);
v    =uniform_flow_v(zv)   +ellipse_flow_v(zv);
// factor potential by alpha coefficients inside circles for continuous equipotential contours across the interface
global ellipse
alpha=ones(z);
for j=1:max(size(ellipse))
  zeta=(z-ellipse(j).zc)/(ellipse(j).Lmajor+ellipse(j).Lminor)*exp(-%i*ellipse(j).alpha);
  Z=zeta+((zeta+ellipse(j).F).^(.5)).*((zeta-ellipse(j).F).^(0.5));
  alpha=alpha+(ellipse(j).solve_alpha_minus/solve_alpha_plus-alpha).*(abs(Z)<1);
end
omega=real(omega).*alpha+%i*imag(omega);

draw_flow_conduction(x,y,real(omega),[0:.5:10],imag(omega),[-100:.5:100],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
ellipse_flow_draw(%F);
