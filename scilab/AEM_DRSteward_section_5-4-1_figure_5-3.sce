clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 5.4.1, Figure 5.3";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

function line_flow_specify(z1,z2,N,M)
// z1,z2: endpoints of line 
// solve_type: type of boundary condition
// N: number of coefficients (0:N); use default if not specified
// M: number of control points where boundary condition is applied; use default if not specified
  global line;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(1.5*(N+1)); end
  if (N>M) then
    str=sprintf("Error: number of coefficients greater than number of unknowns for line %d\n",max(size(line))+1);
    error(str);
  end
  line_ielt.z1=z1;
  line_ielt.z2=z2;
  line_ielt.L=abs(z2-z1)/2;
  line_ielt.expitheta=(z2-z1)/abs(z2-z1);
  line_ielt.coef_double_layer=zeros(N+1,1);
  line_ielt.solve_Z=(2*[1:M]'-M-1)/M;
  line_ielt.solve_z=line_ielt.solve_Z*(line_ielt.z2-line_ielt.z1)/2+(line_ielt.z2+line_ielt.z1)/2;
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

xmin=-10; xmax=10; ymin=-10; ymax=10;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/36:1/18:1]; yv=ymin+(ymax-ymin)*[1/36:1/18:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

z1=-1; z2=1;
line_flow_specify(z1,z2);

// uncomment to display influence functions for series with continuous stream function
global line 
line(1).coef_double_layer(1)=1;
//line(1).coef_double_layer(2)=1;
//line(1).coef_double_layer(3)=1;
//line(1).coef_double_layer(4)=1;

omega=line_flow_omega(z);
v=line_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[-1:1/20:1],imag(omega),[-2:1/20:2],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
line_flow_draw();
