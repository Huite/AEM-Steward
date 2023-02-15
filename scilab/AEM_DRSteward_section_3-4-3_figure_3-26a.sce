clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.4.3, Figure 3.26a";
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
  omega=slit_flow_omega(z);
endfunction

function [change]=uniform_flow_solve()
  global uniform
  omega_beforesolve=uniform_flow_omega(uniform.solve_z);
  if max(size(uniform.solve_z))>2 then
    // compute Phi0 and v0 for more than 2 control points
    A=[ones(uniform.solve_z),-real(uniform.solve_z),-imag(uniform.solve_z)]; 
    b=real(uniform.solve_Phi-omega_add_uniform_flow(uniform.solve_z));
    coef=A\b;
    uniform.coef_Phi0=coef(1);
    uniform.coef_v0=coef(2)+%i*coef(3);
  else
    // compute Phi0 if less than 2 control points
    A=[ones(uniform.solve_z)]; 
    b=real(uniform.solve_Phi-omega_add_uniform_flow(uniform.solve_z));
    coef=A\b;
    uniform.coef_Phi0=coef;
  end
  omega_aftersolve=uniform_flow_omega(uniform.solve_z);
  change=max(abs(real(omega_aftersolve-omega_beforesolve)));
endfunction

function uniform_flow_error()
  global uniform
  z=uniform.solve_z;
  omega=uniform_flow_omega(z)+slit_flow_omega(z);
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

function slit_flow_specify(z1,z2,solve_type,solve_value,N,M)
// z1,z2: endpoints of the slit element
// solve_type: type of boundary condition
// solve_value: boundary condition takes on values associated with
//   solve_type==solve_type_Phi_specified()  potential (enter as a constant and then can change the values at control points later if variation desired)
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
  if sliti.solve_type==solve_type_Phi_specified() then
    sliti.solve_Phi=solve_value;
  end
  sliti.solve_theta=%pi*([1:M]'-0.5)/M;
  sliti.solve_Zp=exp( %i*sliti.solve_theta);
  sliti.solve_Zm=exp(-%i*sliti.solve_theta);
  sliti.solve_zeta=0.5*(sliti.solve_Zp+sliti.solve_Zp.^(-1));
  sliti.solve_z=(sliti.z2+sliti.z1)/2+sliti.solve_zeta*(sliti.z2-sliti.z1)/2;
  // setup matrices to be used in solving system of equations to determine coefficients
  if solve_type==solve_type_Phi_specified() then
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
    if slit(ielt).solve_type==solve_type_Phi_specified() then
      Phi_solve_z=real(omega_add_slit_flow(slit(ielt).solve_z,[ielt]));
      x=2/M*slit(ielt).A'*(slit(ielt).solve_Phi-Phi_solve_z);
      slit(ielt).coef_Q=-x(1)/2;
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
    if slit(ielt).solve_type==solve_type_Phi_specified() then
      slit(ielt).residual_error=real(omega)-slit(ielt).solve_Phi;
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

xmin=-5; xmax=5; ymin=0; ymax=10;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/24:1/12:1]; yv=ymin+(ymax-ymin)*[1/24:1/12:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

solve_z=[xmin+%i*ymin]; solve_Phi=[10];
uniform_flow_specify(solve_z,solve_Phi);

xyznode=[0.0 0.0 0; // node 1
         0.0 1.5 1; // node 2
         0.7 2.2 2; // node 3
        -0.8 2.5 2; // node 4
         0.8 3.2 3; // node 5
         1.9 3.3 4; // node 6
         1.7 4.0 4; // node 7
        -1.1 3.5 3; // node 8
        -2.2 3.8 4; // node 9
        -0.3 4.3 4; // node 10
         0.7 5.2 5; // node 11
         1.7 4.9 6; // node 12
         0.8 6.6 6; // node 13
         2.0 6.5 7; // node 14
         2.5 5.5 8; // node 15
         2.6 7.0 8; // node 16
         1.0 7.8 7; // node 17
         0.6 8.7 8; // node 18
         1.7 8.3 8; // node 19
        -1.2 5.2 5; // node 20
        -2.0 5.0 6; // node 21
        -0.8 6.7 6; // node 22
        -2.0 6.8 7; // node 23
        -2.7 6.1 8; // node 24
        -2.4 7.6 8; // node 25
        -0.8 8.0 7; // node 26
        -1.2 8.7 8; // node 27
        -0.2 8.8 8; // node 28
         ];
endpoints=[1 2;
           2 3;
           2 4;
           4 5;
           5 6;
           5 7;
           4 8;
           8 9;
           8 10;
           10 11;
           11 12;
           11 13;
           13 14;
           14 15;
           14 16;
           13 17;
           17 18;
           17 19;
           10 20;
           20 21;
           20 22;
           22 23;
           23 24;
           23 25;
           22 26;
           26 27;
           26 28;
           ];

solve_type=solve_type_Phi_specified(); 
for ielt=1:size(endpoints,1)
  z1=xyznode(endpoints(ielt,1),1)+%i*xyznode(endpoints(ielt,1),2);
  z2=xyznode(endpoints(ielt,2),1)+%i*xyznode(endpoints(ielt,2),2);
  solve_value=[xyznode(endpoints(ielt,1),3),xyznode(endpoints(ielt,2),3)];
  slit_flow_specify(z1,z2,solve_type,solve_value);
  global slit; slit($).solve_Phi=solve_value(1)*(real(slit($).solve_Zp)-1)/(-2)+solve_value(2)*(real(slit($).solve_Zp)+1)/2;
end

printf("Solving for unknown coefficients\n");
count=0; change=100;
while count<5000 & change>1e-12 do
  count=count+1;
  change=uniform_flow_solve();
  change=max(change,slit_flow_solve());
  printf("Solving iteration %3d change %12.5e\n",count,change);
end
uniform_flow_error();
slit_flow_error();

omega=uniform_flow_omega(z)+slit_flow_omega(z);
v    =uniform_flow_v(zv)   +slit_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[0:0.25:13],[],[],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
uniform_flow_draw();
slit_flow_draw();
