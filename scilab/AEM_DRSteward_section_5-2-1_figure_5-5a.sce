clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 5.2.1, Figure 5.5a";
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

function line_flow_specify(z1,z2,solve_type,solve_value,N,M)
// z1,z2: endpoints of line 
// solve_type: type of boundary condition
// solve_value: boundary condition takes on values associated with
//   solve_type==solve_type_Phi_uniform()    set to [] since unused
// N: number of coefficients (0:N); use default if not specified
// M: number of control points where boundary condition is applied; use default if not specified
  global line;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(1.5*(N+1)); end
  line_ielt.z1=z1;
  line_ielt.z2=z2;
  line_ielt.L=abs(z2-z1)/2;
  line_ielt.expitheta=(z2-z1)/abs(z2-z1);
  line_ielt.coef_double_layer=zeros(N+1,1);
  line_ielt.solve_type=solve_type;
  line_ielt.solve_Z=(2*[1:M]'-M-1)/M;
  line_ielt.solve_z=line_ielt.solve_Z*(line_ielt.z2-line_ielt.z1)/2+(line_ielt.z2+line_ielt.z1)/2;
  if solve_type==solve_type_Phi_uniform() then
    line_ielt.solve_A=ones(line_ielt.solve_Z);
    for n=0:N
      nearfield=zeros(line_ielt.solve_Z);
      for l=1:1:floor(0.5*(n+1))
        nearfield=nearfield-2/(2*l-1)*line_ielt.solve_Z.^(n+1-2*l);
      end
      omega_double=1/(2*%pi)*((line_ielt.solve_Z.^n).*(log(line_ielt.solve_Z+1)-log(line_ielt.solve_Z-1))+nearfield);
      line_ielt.solve_A=[line_ielt.solve_A,real(omega_double)]; 
    end
  end
  if max(size(line))==0 then
    line=list(line_ielt);
  else
    line($+1)=line_ielt;
  end
endfunction

function [omega]=line_flow_omega(z,excluded)
  global line;
  if ~exists('excluded') then excluded=[]; end
  omega=zeros(z);
  for ielt=1:max(size(line))
    if find(ielt==excluded)==[] then 
      // compute local coordinates and remove singularities
      Z=(z-0.5*(line(ielt).z2+line(ielt).z1))./(0.5*(line(ielt).z2-line(ielt).z1));  R=abs(Z);  
      endtol=1e-5; // tolerance for singularities
                       Z=Z+(endtol*%i-Z).*(R<endtol);     // move points near origin to + side of element
      Zend1=-1+endtol; Z=Z+(Zend1-Z).*(abs(Z+1)<endtol);  // move points near Z=-1 onto element
      Zend2= 1-endtol; Z=Z+(Zend2-Z).*(abs(Z-1)<endtol);  // move points near Z=+1 onto element
      // precompute all powers of Z for near and far field
      N=max(size(line(ielt).coef_double_layer))-1;
      Z_power=list(ones(Z)); for n=1:N+1 Z_power($+1)=Z.^(n); end
      // near field double layer
      logZp1=log(Z+1); logZm1=log(Z-1); logZp1overZm1=logZp1-logZm1;
      for n=0:N
        nearfield=zeros(Z);
        for l=1:1:floor(0.5*(n+1))
          nearfield=nearfield-2/(2*l-1)*Z_power(n+1-2*l +1);
        end
        omega=omega+line(ielt).coef_double_layer(n+1)/(2*%pi)*(Z_power(n+1).*logZp1overZm1+nearfield);
      end
    end
  end
endfunction

function [v]=line_flow_v(z,excluded)
   global line;
  if ~exists('excluded') then excluded=[]; end
  v=zeros(z);
  for ielt=1:max(size(line))
    if find(ielt==excluded)==[] then 
      // compute local coordinates and remove singularities
      Z=(z-0.5*(line(ielt).z2+line(ielt).z1))./(0.5*(line(ielt).z2-line(ielt).z1));  R=abs(Z);  
      endtol=1e-5; // tolerance for singularities
                       Z=Z+(endtol*%i-Z).*(R<endtol);     // move points near origin to + side of element
      Zend1=-1+endtol; Z=Z+(Zend1-Z).*(abs(Z+1)<endtol);  // move points near Z=-1 onto element
      Zend2= 1-endtol; Z=Z+(Zend2-Z).*(abs(Z-1)<endtol);  // move points near Z=+1 onto element
      // precompute all powers of Z for near and far field
      N=max(size(line(ielt).coef_double_layer))-1;
      Z_power=list(ones(Z)); for n=1:N+1 Z_power($+1)=Z.^(n); end
      // near field double layer
      V=zeros(Z);
      logZp1=log(Z+1); logZm1=log(Z-1); logZp1overZm1=logZp1-logZm1;
      Zp1_minus1=(Z+1).^(-1); Zm1_minus1=(Z-1).^(-1);
      for n=0:N
        nearfield=zeros(Z);
        for l=1:1:floor(0.5*(n))
          nearfield=nearfield+2/(2*l-1)*Z_power(n-2*l +1);
        end
        Vbefore_conj=((-1)^n)*Zp1_minus1-Zm1_minus1-n*nearfield;  
        if n>0 then
          Vbefore_conj=Vbefore_conj+n*Z_power(n-1 +1).*logZp1overZm1;
        end
        V=V-conj(line(ielt).coef_double_layer(n+1)/(2*%pi)*Vbefore_conj); 
      end
      v=v+V/conj(0.5*(line(ielt).z2-line(ielt).z1));
    end
  end
endfunction

function [omega]=omega_add_line_flow(z,excluded)
  omega=uniform_flow_omega(z)...
       +line_flow_omega(z,excluded);
endfunction

function [change]=line_flow_solve()
  global line
  change=0;
  for ielt=1:max(size(line))
    omega_beforesolve=line_flow_omega(line(ielt).solve_z,[1:ielt-1,ielt+1:size(line)]);
    omega_add=omega_add_line_flow(line(ielt).solve_z,[ielt])
    if line(ielt).solve_type==solve_type_Phi_uniform() then
      line(ielt).solve_b=-real(omega_add);
      coef=lsq(line(ielt).solve_A,line(ielt).solve_b);
      line(ielt).coef_double_layer=coef(2:$);
      line(ielt).solve_Phi=-coef(1);
    end
    omega_aftersolve=line_flow_omega(line(ielt).solve_z,[1:ielt-1,ielt+1:max(size(line))]);
    change=max(change,max(abs(real(omega_aftersolve-omega_beforesolve))));
  end
endfunction

function line_flow_error()
  global line
  if max(size(line))<1 then return end
  for ielt=1:size(line)
    z=line(ielt).solve_z;  Z=line(ielt).solve_Z;
    omega=line_flow_omega(z)+uniform_flow_omega(z);
    v    =line_flow_v(z)    +uniform_flow_v(z)    ;
    if line(ielt).solve_type==solve_type_Phi_uniform() then
      line(ielt).residual_error=real(omega)-line(ielt).solve_Phi;
    end
  end
  net_error=[];
  for ielt=1:size(line)
    error_rmse_ielt = sqrt( sum(    ( line(ielt).residual_error ).^2 )/max(size( line(ielt).residual_error )) );     
    error_abs_ielt  =       sum( abs( line(ielt).residual_error )    )/max(size( line(ielt).residual_error ))  ;
    net_error=[net_error;line(ielt).residual_error];
    printf("Line %3d error rmse: %10.4e abs: %10.4e %s\n",ielt,error_rmse_ielt,error_abs_ielt,line(ielt).solve_type);
  end
  if max(size(net_error))>0 then printf("Line net error rmse: %10.4e average absolute value: %10.4e \n",...
     sqrt( sum(net_error.^2)/max(size(net_error)) ), sum(abs(net_error))/max(size(net_error)) ); end
endfunction

function line_flow_draw(draw_geometry)
  global line;
  if ~exists('draw_geometry') then draw_geometry=%F; end
  drawlater();
  for ielt=1:max(size(line))
    draw_segments(real(line(ielt).z1),imag(line(ielt).z1),real(line(ielt).z2),imag(line(ielt).z2));
    if draw_geometry then
      draw_points(real([line(ielt).z1;line(ielt).z2]),imag([line(ielt).z1;line(ielt).z2]));
      draw_points(real(line(ielt).solve_z),imag(line(ielt).solve_z));
    end
  end
  drawnow();
endfunction

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

// set the uniform flow components for Phi0 and v0 instead of solving to match boundary conditions
uniform_flow_specify([],[]); global uniform; uniform.coef_Phi0=0; uniform.coef_v0=1;

solve_type=solve_type_Phi_uniform(); solve_value=[];
z1=-1-%i; z2=1+%i;
line_flow_specify(z1,z2,solve_type,solve_value);

printf("Solving for unknown coefficients\n");
count=0; change=100;
while count<5000 & change>1e-12 do
  count=count+1;
  change=0;
  change=max(change,line_flow_solve());
  printf("Solving iteration %3d change %12.5e\n",count,change);
end
line_flow_error();

omega=uniform_flow_omega(z)+line_flow_omega(z);
v    =uniform_flow_v(zv)   +line_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[-2:1/10:2],imag(omega),[-10:1/10:10],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
line_flow_draw();
