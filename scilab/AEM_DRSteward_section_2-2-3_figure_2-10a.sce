clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 2.2.1, Figure 2.10a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

s_m=[0:1:4]'; s_plot=[min(s_m):(max(s_m)-min(s_m))/100:max(s_m)]'; M=max(size(s_m));
smin=min(s_m); smax=max(s_m); 
fmin=0; fmax=2;
f_m=[0.04;0.38;1.03;1.77;1.96];
c=[0;0]; csave=[[-100;-100],c];
while max(abs(csave(:,$)-csave(:,$-1)))>1e-12 
  fapprox=fmin+(fmax-fmin)*(1+exp(c(1)+c(2)*s_m)).^(-1);
  J=[-(fmax-fmin)*   exp(c(1)+c(2)*s_m).*(1+exp(c(1)+c(2)*s_m)).^(-2), ...
     -(fmax-fmin)*s_m.*exp(c(1)+c(2)*s_m).*(1+exp(c(1)+c(2)*s_m)).^(-2)];
  cnext=c-J\(fapprox-f_m);
  csave=[csave,cnext]; c=cnext;
end
printf("c0:%f c1:%f\n",c(1),c(2));
f_approx=fmin+(fmax-fmin)*(1+exp(c(1)+c(2)*s_m)).^(-1);
printf("objective:%f\n",sum((f_approx-f_m).^2)/M);
f_plot=fmin+(fmax-fmin)*(1+exp(c(1)+c(2)*s_plot)).^(-1);
draw_function(s_plot,f_plot,s_m,f_m);
xtitle(title_author,"x","f");
