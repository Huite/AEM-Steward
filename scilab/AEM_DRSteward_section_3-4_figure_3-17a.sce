clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.4, Figure 3.17a";
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
  if max(size(slit))==0 then
    slit=list(sliti);
  else
    slit($+1)=sliti;
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

xmin=0; xmax=100; ymin=0; ymax=100;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

z1=30+%i*15; z2=70+%i*85;
slit_flow_specify(z1,z2);

draw_flow_conduction(x,y,[],[],[],[],xv,yv,[],[]);
xtitle(title_author,"x","y");
slit_flow_draw(%T);
