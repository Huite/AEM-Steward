clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 5.2.1, Figure 5.6c";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// type of boundary/interface condition
function [solve_type]=solve_type_Robin_dPsi()
  solve_type="solution type: Robin tangential vector from discontinuous stream function"; 
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
//   solve_type==solve_type_Robin_dPsi() delta = line k * line width / regional k
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
  if solve_type==solve_type_Robin_dPsi() then
    line_ielt.solve_delta=solve_value;
    line_ielt.solve_A=[];
    for n=0:N
      nearfield=zeros(line_ielt.solve_Z);
      for l=1:1:floor(0.5*(n))
        nearfield=nearfield-2/(2*l-1)*line_ielt.solve_Z.^(n-2*l);
      end
      vsn_double=-1/(2*%pi*line_ielt.L)*conj((-1)^n*(line_ielt.solve_Z+1).^(-1)-(line_ielt.solve_Z-1).^(-1)+n*nearfield);
      if n>0 then
        vsn_double=vsn_double-1/(2*%pi*line_ielt.L)*conj(n*(line_ielt.solve_Z.^(n-1)).*(log(line_ielt.solve_Z+1)-log(line_ielt.solve_Z-1)));
      end
      line_ielt.solve_A=[line_ielt.solve_A,real(vsn_double)-(real(line_ielt.solve_Z)).^n/line_ielt.solve_delta]; 
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
  R_far=1.5; N_far=12; // far field tolerance and number of terms for far field expansions
  for ielt=1:max(size(line))
    if find(ielt==excluded)==[] then 
      // compute local coordinates and remove singularities
      Z=(z-0.5*(line(ielt).z2+line(ielt).z1))./(0.5*(line(ielt).z2-line(ielt).z1));  R=abs(Z);  
      endtol=1e-5; // tolerance for singularities
                       Z=Z+(endtol*%i-Z).*(R<endtol);     // move points near origin to + side of element
      Zend1=-1+endtol; Z=Z+(Zend1-Z).*(abs(Z+1)<endtol);  // move points near Z=-1 onto element
      Zend2= 1-endtol; Z=Z+(Zend2-Z).*(abs(Z-1)<endtol);  // move points near Z=+1 onto element
      // setup sparse matrices to only compute in either nearfield or farfield
      [near_ij,Z_near,near_mn]=spget(sparse(Z.*(R< R_far)));
      [far_ij, Z_far, far_mn ]=spget(sparse(Z.*(R>=R_far)));
      // precompute all powers of Z for near and far field
      N_double=max(size(line(ielt).coef_double_layer))-1;
      N=N_double;
      Zpow_near=list(ones(Z_near)); for n=1:N+1              Zpow_near($+1)=Z_near.^(n); end
      Zpow_far =list(Z_far.^(-1));  for n=-2:-1:-2*(N_far+1) Zpow_far($+1) =Z_far.^(n); end
      if near_ij~=[] then
        // near field double layer
        omega_near=zeros(Z_near);
        logZp1=log(Z_near+1);  logZm1=log(Z_near-1);  logZp1overZm1=logZp1-logZm1;
        for n=0:N_double
          nearfield=zeros(Z_near);
          for l=1:1:floor(0.5*(n+1))
            nearfield=nearfield-2/(2*l-1)*Zpow_near(n+1-2*l +1);
          end
          omega_near=omega_near+line(ielt).coef_double_layer(n+1)/(2*%pi)*(Zpow_near(n+1).*logZp1overZm1+nearfield);
        end
      end
      if far_ij~=[] then
        // far field double layer
        omega_far=zeros(Z_far);
        for n=0:N_double
          farfield=zeros(Z_far);
          lrange=floor(0.5*(n+3))+[0:N_far-1];
          for l=lrange
            farfield=farfield+2/(2*l-1)*Zpow_far(-(n+1-2*l));  //farfield=farfield+2/(2*l-1)*Zfar.^(n+1-2*l);
          end
          omega_far=omega_far+line(ielt).coef_double_layer(n+1)/(2*%pi)*farfield;
        end
      end
      if near_ij~=[] then omega=omega+full(sparse(near_ij,omega_near,near_mn)); end
      if far_ij~=[]  then omega=omega+full(sparse(far_ij ,omega_far ,far_mn )); end
    end
  end
endfunction

function [v]=line_flow_v(z,excluded)
   global line;
  if ~exists('excluded') then excluded=[]; end
  v=zeros(z);
  R_far=1.5; N_far=12; // far field tolerance and number of terms for far field expansions
  for ielt=1:max(size(line))
    if find(ielt==excluded)==[] then 
      // compute local coordinates and remove singularities
      Z=(z-0.5*(line(ielt).z2+line(ielt).z1))./(0.5*(line(ielt).z2-line(ielt).z1));  R=abs(Z);  
      endtol=1e-5; // tolerance for singularities
                       Z=Z+(endtol*%i-Z).*(R<endtol);     // move points near origin to + side of element
      Zend1=-1+endtol; Z=Z+(Zend1-Z).*(abs(Z+1)<endtol);  // move points near Z=-1 onto element
      Zend2= 1-endtol; Z=Z+(Zend2-Z).*(abs(Z-1)<endtol);  // move points near Z=+1 onto element
      // setup sparse matrices to only compute in either nearfield or farfield
      [near_ij,Z_near,near_mn]=spget(sparse(Z.*(R< R_far)));
      [far_ij, Z_far, far_mn ]=spget(sparse(Z.*(R>=R_far)));
      // precompute all powers of Z for near and far field
      N_double=max(size(line(ielt).coef_double_layer))-1;
      N=N_double;
      Zpow_near=list(ones(Z_near)); for n=1:N+1              Zpow_near($+1)=Z_near.^(n); end
      Zpow_far =list(Z_far.^(-1));  for n=-2:-1:-2*(N_far+1) Zpow_far($+1) =Z_far.^(n); end
      if near_ij~=[] then
        // near field double layer
        V_near=zeros(Z_near);
        logZp1=log(Z_near+1);  logZm1=log(Z_near-1);  logZp1overZm1=logZp1-logZm1;
        Zp1_minus1=(Z_near+1).^(-1); Zm1_minus1=(Z_near-1).^(-1);
        for n=0:N_double
          nearfield=zeros(Z_near);
          for l=1:1:floor(0.5*(n))
            nearfield=nearfield+2/(2*l-1)*Zpow_near(n-2*l +1);
          end
          Vbefore_conj=((-1)^n)*Zp1_minus1-Zm1_minus1-n*nearfield;  
          if n>0 then
            Vbefore_conj=Vbefore_conj+n*Zpow_near(n-1 +1).*logZp1overZm1;
          end
          V_near=V_near-conj(line(ielt).coef_double_layer(n+1)/(2*%pi)*Vbefore_conj); 
        end
      end
      if far_ij~=[] then
        // far field double layer
        V_far=zeros(Z_far);
        for n=0:N_double
          farfield=zeros(Z_far);
          lrange=floor(0.5*(n+3))+[0:N_far-1];
          for l=lrange
            farfield=farfield+2/(2*l-1)*(n+1-2*l)*Zpow_far(-(n-2*l));  
          end
          Vbefore_conj=farfield;  
          V_far=V_far-conj(line(ielt).coef_double_layer(n+1)/(2*%pi)*Vbefore_conj); 
        end
      end
      V=0;
      if near_ij~=[] then V=V+full(sparse(near_ij,V_near,near_mn)); end
      if far_ij~=[]  then V=V+full(sparse(far_ij ,V_far ,far_mn ));; end
      v=v+V/conj(0.5*(line(ielt).z2-line(ielt).z1));
    end
  end
endfunction

// additional function without the excluded line element
function [v]=v_add_line_flow(z,excluded)
  v=uniform_flow_v(z)...
   +line_flow_v(z,excluded);
endfunction

function [change]=line_flow_solve()
  global line
  change=0;
  for ielt=1:max(size(line))
    omega_beforesolve=line_flow_omega(line(ielt).solve_z,[1:ielt-1,ielt+1:size(line)]);
    v_add=v_add_line_flow(line(ielt).solve_z,[ielt]);
    if line(ielt).solve_type==solve_type_Robin_dPsi() then
      line(ielt).solve_b=-real(v_add*conj(line(ielt).expitheta));
      coef=lsq(line(ielt).solve_A,line(ielt).solve_b);
      line(ielt).coef_double_layer=coef;
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
    if line(ielt).solve_type==solve_type_Robin_dPsi() then
      delta_Psi=0;
      for n=0:max(size(line(ielt).coef_double_layer))-1
        delta_Psi=delta_Psi-real(line(ielt).coef_double_layer(n+1))*(real(Z)).^n;
      end
      line(ielt).residual_error=real(v*conj(line(ielt).expitheta))+delta_Psi/line(ielt).solve_delta;
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

xmin=-5; xmax=5; ymin=0; ymax=10;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

solve_z=[xmin+%i*ymax;xmin+%i*ymin;xmax+%i*ymin;xmax+%i*ymax]; solve_Phi=[1;1;-1;-1];
uniform_flow_specify(solve_z,solve_Phi);

// compute random locations along non-intersecting adjacent families of line elements
J=4; Lmin=(ymax-ymin)/5; Lmax=(ymax-ymin)/3;
z1save=[]; z2save=[];
for count=1:J
  overlap=%T;
  while overlap do 
    istring=3;
    zstart=grand(1,1,'unf',xmin,xmin+Lmin)+%i*grand(1,1,'unf',ymin+Lmin,ymax-Lmin);
    L=grand(istring,1,'unf',Lmin,Lmax);
    angle=grand(istring,1,'unf',-%pi/4,%pi/4);
    z1test=zstart; z2test=z1test+L(1)*exp(%i*angle(1));
    for k=2:istring
      z1test(k)=z2test(k-1); z2test(k)=z1test(k)+L(k)*exp(%i*angle(k));
    end
    overlap=%F;
    for j=1:max(size(z1save))
      for k=1:istring
        z1=z1test(k); z2=z2test(k);
        Zend=([z1save(j);z2save(j)]-(z2+z1)/2)*2/(z2-z1);
        X1=real(Zend(1)); X2=real(Zend(2)); Y1=imag(Zend(1)); Y2=imag(Zend(2));
        tol=0.02*(ymax-ymin); // tol=0.05;
        if -1-tol<X1 & X1<1+tol & abs(Y1)<tol then overlap=%T; end
        if -1-tol<X2 & X2<1+tol & abs(Y2)<tol then overlap=%T; end
        if Y1*Y2<=0 then
          Xintercept=X2-Y2*(X2-X1)/(Y2-Y1);
          if -1-tol<Xintercept & Xintercept<1+tol then
            overlap=%T;
          end
        end
        Zend=([z1;z2]-(z2save(j)+z1save(j))/2)*2/(z2save(j)-z1save(j));
        X1=real(Zend(1)); X2=real(Zend(2)); Y1=imag(Zend(1)); Y2=imag(Zend(2));
        if -1-tol<X1 & X1<1+tol & abs(Y1)<tol then overlap=%T; end
        if -1-tol<X2 & X2<1+tol & abs(Y2)<tol then overlap=%T; end
      end
    end
    if min(imag(z2test))<ymin | max(imag(z2test))>ymax then 
      overlap=%T;
    end
    if ~overlap then
      z1save=[z1save;z1test]; z2save=[z2save;z2test];
    end
  end
end
for i=1:max(size(z1save))
  solve_type=solve_type_Robin_dPsi(); solve_value=1;
  line_flow_specify(z1save(i),z2save(i),solve_type,solve_value);
end

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
line_flow_draw();
