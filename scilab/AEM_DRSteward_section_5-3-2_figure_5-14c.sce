clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 5.3.2, Figure 5.14c";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

// type of boundary/interface condition
function [solve_type]=solve_type_Kutta()
  solve_type="solution type: Kutta condition and stream function uniform";
endfunction
function [solve_type]=solve_type_Psi_specified()
  solve_type="solution type: stream function specified"; 
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
//   solve_type==solve_type_Kutta()          %T: left endpoint, %F: right endpoint (default %T if [])
//   solve_type==solve_type_Psi_specified()  stream function
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
  line_ielt.compute_add_single_omega=[];
  if solve_type==solve_type_Kutta() then
    if solve_value==[] then solve_value=%T; end
    line_ielt.solve_Kutta_end1=solve_value;
    line_ielt.solve_Psi=0;
    line_ielt.solve_A=[];
    for n=0:N
      nearfield=zeros(line_ielt.solve_Z);
      for l=1:1:floor(0.5*(n+2))
        nearfield=nearfield-2/(2*l-1)*line_ielt.solve_Z.^(n+2-2*l);
      end
      omega_single=1/(2*%pi*(n+1))*((line_ielt.solve_Z.^(n+1)).*(log(line_ielt.solve_Z+1)-log(line_ielt.solve_Z-1))...
                                                   +((-1)^n)*log(line_ielt.solve_Z+1)+log(line_ielt.solve_Z-1)+nearfield);
      line_ielt.solve_A=[line_ielt.solve_A,imag(%i*omega_single)]; 
    end
    if line_ielt.solve_Kutta_end1 then
      line_ielt.solve_A=[line_ielt.solve_A;(-1).^[0:N]];
    else
      line_ielt.solve_A=[line_ielt.solve_A;(+1).^[0:N]];
    end
    line_ielt.solve_A=[[ones(M,1);0],line_ielt.solve_A];
  elseif solve_type==solve_type_Psi_specified() then
    line_ielt.solve_Psi=solve_value;  // constant value of Phi along element
    line_ielt.solve_A=[];
    for n=0:N
      nearfield=zeros(line_ielt.solve_Z);
      for l=1:1:floor(0.5*(n+2))
        nearfield=nearfield-2/(2*l-1)*line_ielt.solve_Z.^(n+2-2*l);
      end
      omega_single=1/(2*%pi*(n+1))*((line_ielt.solve_Z.^(n+1)).*(log(line_ielt.solve_Z+1)-log(line_ielt.solve_Z-1))...
                                                   +((-1)^n)*log(line_ielt.solve_Z+1)+log(line_ielt.solve_Z-1)+nearfield);
      line_ielt.solve_A=[line_ielt.solve_A,imag(%i*omega_single)]; 
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
      N_single=max(size(line(ielt).coef_single_layer))-1;
      N=N_single;
      Zpow_near=list(ones(Z_near)); for n=1:N+1              Zpow_near($+1)=Z_near.^(n); end
      Zpow_far =list(Z_far.^(-1));  for n=-2:-1:-2*(N_far+1) Zpow_far($+1) =Z_far.^(n); end
      if near_ij~=[] then 
        // near field single layer
        omega_near=zeros(Z_near);
        logZp1=log(Z_near+1);  logZm1=log(Z_near-1);  logZp1overZm1=logZp1-logZm1;
        for n=0:N_single
          nearfield=zeros(Z_near);
          for l=1:1:floor(0.5*(n+2))
            nearfield=nearfield-2/(2*l-1)*Zpow_near(n+2-2*l +1);
          end
          omega_near=omega_near+line(ielt).coef_single_layer(n+1)/(2*%pi)/(n+1)...
                     *(Zpow_near(n+2).*logZp1overZm1+(-1)^n*logZp1+logZm1+nearfield);
        end
      end
      if far_ij~=[] then 
        // far field single layer
        omega_far=zeros(Z_far);
        logZp1_far=log(Z_far+1); logZm1_far=log(Z_far-1);
        for n=0:N_single
          farfield=zeros(Z_far);
          lrange=floor(0.5*(n+4))+[0:N_far-1];
          for l=lrange
            farfield=farfield+2/(2*l-1)*Zpow_far(-(n+2-2*l));  //farfield=farfield+2/(2*l-1)*Zfar.^(n+2-2*l);
          end
          omega_far=omega_far+line(ielt).coef_single_layer(n+1)/(2*%pi)/(n+1)*((-1)^n*logZp1_far+logZm1_far+farfield);
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
      N_single=max(size(line(ielt).coef_single_layer))-1;
      N=N_single;
      Zpow_near=list(ones(Z_near)); for n=1:N+1              Zpow_near($+1)=Z_near.^(n); end
      Zpow_far =list(Z_far.^(-1));  for n=-2:-1:-2*(N_far+1) Zpow_far($+1) =Z_far.^(n); end
      if near_ij~=[] then 
         // near field single layer
        V_near=zeros(Z_near);
        logZp1=log(Z_near+1);  logZm1=log(Z_near-1);  logZp1overZm1=logZp1-logZm1;
        for n=0:N_single
          nearfield=zeros(Z_near);
          for l=1:floor(0.5*(n+1))
            nearfield=nearfield-2/(2*l-1)*Zpow_near(n+1-2*l +1);
          end
          Vbefore_conj=Zpow_near(n +1).*logZp1overZm1+nearfield;  
          V_near=V_near-conj(line(ielt).coef_single_layer(n+1)/(2*%pi)*Vbefore_conj); 
        end
      end
      if far_ij~=[] then 
        // far field single layer
        V_far=zeros(Z_far);
        for n=0:N_single
          farfield=zeros(Z_far);
          lrange=floor(0.5*(n+3))+[0:N_far-1];
          for l=lrange
            farfield=farfield-2/(2*l-1)*Zpow_far(-(n+1-2*l));  //farfield=farfield-2/(2*l-1)*Zfar.^(n+1-2*l);
          end
          Vbefore_conj=farfield;  
          V_far=V_far+conj(line(ielt).coef_single_layer(n+1)/(2*%pi)*Vbefore_conj); 
        end
      end
      V=0;
      if near_ij~=[] then V=V+full(sparse(near_ij,V_near,near_mn)); end
      if far_ij~=[]  then V=V+full(sparse(far_ij ,V_far ,far_mn ));; end
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
  // fill matrices used to compute additional functions with matrix multiplication instead of function evaluation
  compute_add_matrices=%T;
  if compute_add_matrices then
    if line(1).compute_add_single_omega==[] then
      printf("Precomputing additional element matrix multiplication ");
      for ielt=1:max(size(line))
        printf(".");
        solve_z=line(ielt).solve_z;
        for jelt=1:size(line)
          compute_add_single_omega=[];
          coef_single_layer=line(jelt).coef_single_layer; line(jelt).coef_single_layer=zeros(line(jelt).coef_single_layer);
          for n=0:max(size(coef_single_layer))-1
            line(jelt).coef_single_layer(n+1)=1;
            compute_add_single_omega=[compute_add_single_omega,line_flow_omega(solve_z,[1:jelt-1,jelt+1:size(line)])];
            line(jelt).coef_single_layer(n+1)=0;
          end
          line(jelt).coef_single_layer=coef_single_layer;
          if jelt==1 then 
            line(ielt).compute_add_single_omega=list(compute_add_single_omega); 
          else 
            line(ielt).compute_add_single_omega($+1)=compute_add_single_omega; 
          end
        end
      end
      printf("\n");
    end
  end
  // Solve coefficients for each line element 
  for ielt=1:max(size(line))
    omega_beforesolve=line_flow_omega(line(ielt).solve_z,[1:ielt-1,ielt+1:size(line)]);
    // compute additional functions
    if compute_add_matrices then
      omega_add=zeros(line(ielt).solve_z); 
      jelt_range=[1:ielt-1,ielt+1:size(line)];
      for jelt=jelt_range
        omega_add=omega_add+line(ielt).compute_add_single_omega(jelt)*line(jelt).coef_single_layer;
      end
      omega_add=omega_add+uniform_flow_omega(line(ielt).solve_z);
    else
      omega_add=omega_add_line_flow(line(ielt).solve_z,[ielt])
    end
    if line(ielt).solve_type==solve_type_Kutta() then
      line(ielt).solve_b=[-imag(omega_add);0];
      coef=lsq(line(ielt).solve_A,line(ielt).solve_b);
      line(ielt).solve_Psi=-coef(1);
      line(ielt).coef_single_layer=%i*coef(2:$);
    elseif line(ielt).solve_type==solve_type_Psi_specified() then
      line(ielt).solve_b=line(ielt).solve_Psi-imag(omega_add);
      coef=lsq(line(ielt).solve_A,line(ielt).solve_b);
      SOR=.01; // successive overrelaxation coefficient, fraction of change to apply to each iteration
      line(ielt).coef_single_layer=line(ielt).coef_single_layer+SOR*(%i*coef-line(ielt).coef_single_layer);
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
    if line(ielt).solve_type==solve_type_Kutta() then
      line(ielt).residual_error=imag(omega)-line(ielt).solve_Psi;
    elseif line(ielt).solve_type==solve_type_Psi_specified() then
      line(ielt).residual_error=imag(omega)-line(ielt).solve_Psi;
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

uniform_flow_specify([],[]); global uniform; uniform.coef_Phi0=0; uniform.coef_v0=1*exp(%i*%pi);

L=1; alpha=9.5/180*%pi; N=10;
zwing=L*(-%i/tan(alpha)+1/sin(alpha)*exp(%i*(%pi/2+[alpha:-alpha*2/N:-alpha]')));
zwing1=zwing-%i*.5;
iline=1; line_flow_specify(zwing1(iline),zwing1(iline+1),solve_type_Kutta(),%T);
for iline=2:max(size(zwing1))-1
  line_flow_specify(zwing1(iline),zwing1(iline+1),solve_type_Psi_specified(),3);
end
zwing2=zwing+%i*.5;
iline=1; line_flow_specify(zwing2(iline),zwing2(iline+1),solve_type_Kutta(),%T);
for iline=2:max(size(zwing2))-1
  line_flow_specify(zwing2(iline),zwing2(iline+1),solve_type_Psi_specified(),3);
end

printf("Solving for unknown coefficients\n");
count=0; change=100;
while count<5000 & change>1e-12 do
  count=count+1;
  change=0;
  global line
  for jelt=2:size(line)/2
    line(jelt).solve_Psi=line(1).solve_Psi;
  end
  for jelt=size(line)/2+2:size(line)
    line(jelt).solve_Psi=line(size(line)/2+1).solve_Psi;
  end
  change=line_flow_solve();
  printf("Solving iteration %3d change %12.5e\n",count,change);
end
line_flow_error();

omega=uniform_flow_omega(z)+line_flow_omega(z);
v    =uniform_flow_v(zv)   +line_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[-2:1/10:2],imag(omega),[-10:1/10:10],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
line_flow_draw();
