clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 5.3.1, Figure 5.12a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// type of boundary/interface condition
function [solve_type]=solve_type_Phi_resistance()
  solve_type="solution type: potential specified resistance"; 
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

// additional function without the uniform flow
function [omega]=omega_add_uniform_flow(z)
  omega=line_flow_omega(z);
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
  omega=uniform_flow_omega(z)+line_flow_omega(z);
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

function line_flow_specify(z1,z2,solve_type,solve_value,N,M)
// z1,z2: endpoints of line 
// solve_type: type of boundary condition
// solve_value: boundary condition takes on values associated with
//   solve_type==solve_type_Phi_resistance() resistance, potential (enter as a constant or value at two endpoints)
// N: number of coefficients (0:N); use default if not specified
// M: number of control points where boundary condition is applied; use default if not specified
  global line;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(1.5*(N+1)); end
  line_ielt.z1=z1;
  line_ielt.z2=z2;
  line_ielt.L=abs(z2-z1)/2;
  line_ielt.expitheta=(z2-z1)/abs(z2-z1);
  line_ielt.coef_single_layer=zeros(N+1,1); 
  line_ielt.solve_type=solve_type;
  line_ielt.solve_Z=(2*[1:M]'-M-1)/M;
  line_ielt.solve_z=line_ielt.solve_Z*(line_ielt.z2-line_ielt.z1)/2+(line_ielt.z2+line_ielt.z1)/2;
  if solve_type==solve_type_Phi_resistance() then
    line_ielt.solve_delta=solve_value(1);
    if max(size(solve_value))==2 then
      line_ielt.solve_Phi=solve_value(2)*ones(line_ielt.solve_Z);  // constant value of Phi along element
    else
      line_ielt.solve_Phi_end1=solve_value(2); line_ielt.solve_Phi_end2=solve_value(3);  // linear variation between two values of Phi at each end
      line_ielt.solve_Phi=line_ielt.solve_Phi_end1*(line_ielt.solve_Z-1)/(-2)+line_ielt.solve_Phi_end2*(line_ielt.solve_Z+1)/2;
    end
    line_ielt.solve_A=[];
    for n=0:N
      nearfield=zeros(line_ielt.solve_Z);
      for l=1:1:floor(0.5*(n+2))
        nearfield=nearfield-2/(2*l-1)*line_ielt.solve_Z.^(n+2-2*l);
      end
      omega_single=1/(2*%pi*(n+1))*((line_ielt.solve_Z.^(n+1)).*(log(line_ielt.solve_Z+1)-log(line_ielt.solve_Z-1))...
                                                   +((-1)^n)*log(line_ielt.solve_Z+1)+log(line_ielt.solve_Z-1)+nearfield);
      line_ielt.solve_A=[line_ielt.solve_A,real(omega_single)-(real(line_ielt.solve_Z)).^n*line_ielt.solve_delta/line_ielt.L]; 
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
      N=max(size(line(ielt).coef_single_layer))-1;
      Z_power=list(ones(Z)); for n=1:N+1 Z_power($+1)=Z.^(n); end
      // near field single layer
      logZp1=log(Z+1); logZm1=log(Z-1); logZp1overZm1=logZp1-logZm1;
      for n=0:N
        nearfield=zeros(Z);
        for l=1:1:floor(0.5*(n+2))
          nearfield=nearfield-2/(2*l-1)*Z_power(n+2-2*l +1);
        end
        omega=omega+line(ielt).coef_single_layer(n+1)/(2*%pi)/(n+1)...
                   *(Z_power(n+2).*logZp1overZm1+(-1)^n*logZp1+logZm1+nearfield);
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
      N=max(size(line(ielt).coef_single_layer))-1;
      Z_power=list(ones(Z)); for n=1:N+1 Z_power($+1)=Z.^(n); end
      // near field single layer
      V=zeros(Z);
      logZp1=log(Z+1); logZm1=log(Z-1); logZp1overZm1=logZp1-logZm1;
      for n=0:N
        nearfield=zeros(Z);
        for l=1:floor(0.5*(n+1))
          nearfield=nearfield-2/(2*l-1)*Z_power(n+1-2*l +1);
        end
        Vbefore_conj=Z_power(n +1).*logZp1overZm1+nearfield;  
        V=V-conj(line(ielt).coef_single_layer(n+1)/(2*%pi)*Vbefore_conj); 
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
    if line(ielt).solve_type==solve_type_Phi_resistance() then
      line(ielt).solve_b=line(ielt).solve_Phi-real(omega_add);
      coef=lsq(line(ielt).solve_A,line(ielt).solve_b);
      line(ielt).coef_single_layer=coef;
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
    if line(ielt).solve_type==solve_type_Phi_resistance() then
      delta_vn=0;
      for n=0:max(size(line(ielt).coef_single_layer))-1
        delta_vn=delta_vn-real(line(ielt).coef_single_layer(n+1))*(real(Z)).^n/line(ielt).L;
      end
      line(ielt).residual_error=line(ielt).solve_delta*delta_vn-(line(ielt).solve_Phi-real(omega));
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

uniform_flow_specify([xmin+%i*ymax],[0.5]);

solve_type=solve_type_Phi_resistance(); solve_value=[1,-0.5];
z1=-1-%i; z2=1+%i;
line_flow_specify(z1,z2,solve_type,solve_value);

printf("Solving for unknown coefficients\n");
count=0; change=100;
while count<5000 & change>1e-12 do
  count=count+1;
  change=uniform_flow_solve();
  change=max(change,line_flow_solve());
  printf("Solving iteration %3d change %12.5e\n",count,change);
end
uniform_flow_error();
line_flow_error();

omega=uniform_flow_omega(z)+line_flow_omega(z);
v    =uniform_flow_v(zv)   +line_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[-2:1/10:2],imag(omega),[-10:1/10:10],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
uniform_flow_draw();
line_flow_draw();
