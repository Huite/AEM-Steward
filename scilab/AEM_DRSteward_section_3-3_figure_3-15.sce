clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.3, Figure 3.15";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

function ellipse_flow_specify(zc,Lmajor,Lminor,angle,N,M)
// zc: center of ellipse
// Lmajor,Lminor: 1/2 length of major and minor axes of ellipse
// angle: orientation of ellipse with respect to x-axis in units of radians
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

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/36:1/18:1]; yv=ymin+(ymax-ymin)*[1/36:1/18:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

zc=0; Lmajor=1; Lminor=0.6; angle=0;
ellipse_flow_specify(zc,Lmajor,Lminor,angle);

// uncomment to display influence functions for series with continuous stream function
global ellipse
//ellipse(1).coef_in(1)=1;  ellipse(1).coef_out(1)=0; xv=0; yv=0; zv=0; v=zeros(zv);
ellipse(1).coef_in(2)=-1/(1-ellipse(1).F^(2*1));  ellipse(1).coef_out(1)=1;
//ellipse(1).coef_in(2)=%i/(1+ellipse(1).F^(2*1));  ellipse(1).coef_out(1)=%i;
//ellipse(1).coef_in(3)=-1/(1-ellipse(1).F^(2*2));  ellipse(1).coef_out(2)=1;
//ellipse(1).coef_in(3)=%i/(1+ellipse(1).F^(2*2));  ellipse(1).coef_out(2)=%i;
//ellipse(1).coef_in(4)=-1/(1-ellipse(1).F^(2*3));  ellipse(1).coef_out(3)=1;
//ellipse(1).coef_in(4)=%i/(1+ellipse(1).F^(2*3));  ellipse(1).coef_out(3)=%i;

omega=ellipse_flow_omega(z);
v    =ellipse_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[-2:.2:2],imag(omega),[-2:.2:2],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
ellipse_flow_draw(%F);
