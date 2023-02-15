clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.4.2, Figure 3.23";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

function slit_flow_specify(z1,z2,N,M)
// z1,z2: endpoints of the slit element
// N: number of strength coefficients
// M: number of control points where boundary conditions are applied
  global slit;
  if ~exists('N') then N=20; end
  if ~exists('M') then M=ceil(3*N/2); end
  sliti.z1=z1;
  sliti.z2=z2;
  sliti.L=abs(z2-z1)/2;
  sliti.expminusialpha=conj(z2-z1)/abs(z2-z1);
  sliti.solve_theta=%pi*([1:M]'-0.5)/M;
  sliti.solve_Zp=exp( %i*sliti.solve_theta);
  sliti.solve_Zm=exp(-%i*sliti.solve_theta);
  sliti.solve_zeta=0.5*(sliti.solve_Zp+sliti.solve_Zp.^(-1));
  sliti.solve_z=(sliti.z2+sliti.z1)/2+sliti.solve_zeta*(sliti.z2-sliti.z1)/2;
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

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/36:1/18:1]; yv=ymin+(ymax-ymin)*[1/36:1/18:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

z1=-1; z2=1;
slit_flow_specify(z1,z2);

// uncomment to display influence functions for series with continuous stream function
global slit
slit(1).coef(1)=%i;
//slit(1).coef(2)=%i;
//slit(1).coef(3)=%i;
//slit(1).coef(4)=%i;

omega=slit_flow_omega(z);
v    =slit_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[-1:.1:1],imag(omega),[-2:.1:2],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
slit_flow_draw();
