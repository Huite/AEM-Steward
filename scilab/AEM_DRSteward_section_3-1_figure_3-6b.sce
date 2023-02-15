clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 3.1, Figure 3.6b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

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

function []=uniform_flow_solve()
  global uniform
  A=[ones(uniform.solve_z),-real(uniform.solve_z),-imag(uniform.solve_z)]; 
  b=real(uniform.solve_Phi);
  coef=A\b;
  uniform.coef_Phi0=coef(1);
  uniform.coef_v0=coef(2)+%i*coef(3);
endfunction

function uniform_flow_error()
  global uniform
  z=uniform.solve_z;
  omega=uniform_flow_omega(z);
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

xmin=0; xmax=100; ymin=0; ymax=100;
x =xmin+(xmax-xmin)*[0:1/401:1];   y= ymin+(ymax-ymin)*[0:1/401:1]';   z=(x'*ones(y'))'+%i*(ones(x')*y')';
xv=xmin+(xmax-xmin)*[1/20:1/10:1]; yv=ymin+(ymax-ymin)*[1/20:1/10:1]'; zv=(xv'*ones(yv'))'+%i*(ones(xv')*yv')';

solve_z=[10+%i*20;40+%i*80;90+%i*60];
solve_Phi=[20;14;10];
uniform_flow_specify(solve_z,solve_Phi);
uniform_flow_solve();
uniform_flow_error();

omega=uniform_flow_omega(z);
v    =uniform_flow_v(zv);

draw_flow_conduction(x,y,real(omega),[7:1:22],imag(omega),[-300:1:300],xv,yv,real(v),imag(v));
xtitle(title_author,"x","y");
uniform_flow_draw();
