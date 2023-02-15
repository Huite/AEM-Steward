clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 5.1, Figure 5.2a";
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
  line_ielt.solve_Z=(2*[1:M]'-M-1)/M;
  line_ielt.solve_z=line_ielt.solve_Z*(line_ielt.z2-line_ielt.z1)/2+(line_ielt.z2+line_ielt.z1)/2;
  if max(size(line))==0 then
    line=list(line_ielt);
  else
    line($+1)=line_ielt;
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

xmin=-1.5; xmax=1.5; ymin=-1.5; ymax=1.5;
x =xmin+(xmax-xmin)*[0:1/400:1];   y= ymin+(ymax-ymin)*[0:1/400:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

z1=-1.2+%i*0.6;
z2=-0.5+%i*0.9;
z3= 0.5+%i*0.85;
z4= 1.2+%i*0.6;
z5=-1.2-%i*1.1;
z6= 0.5-%i*0.2;
line_flow_specify(z1,z2,5,8);
line_flow_specify(z2,z3,5,8);
line_flow_specify(z3,z4,5,8);
line_flow_specify(z5,z6);

draw_flow_conduction(x,y,[],[],[],[],xv,yv,[],[]);
xtitle(title_author,"x","y");
line_flow_draw(%T);
