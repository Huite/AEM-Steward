clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 2.2.1, Equation (2.34)";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

s_m=[0:1:10]'; s_plot=[0:.1:10]';
f_m=[1.58;1.72;1.70;1.60;0.83;0.46;0.14;0.30;0.51;0.88;1.28]; 
s_exact=10; f_exact=1; 
printf("exact constraint: f(%f)=%f\n",s_exact,f_exact);
M=max(size(s_m));
N=3;
A=ones(s_m); for n=1:N A=[A,s_m.^n]; end
A_exact=ones(s_exact); for n=1:N A_exact=[A_exact,s_exact.^n]; end
coef=[A'*A,A_exact';A_exact,0]\[A'*f_m;f_exact];
c=coef(1:$-1); lambda=coef($);
printf("N:%d\n",N); for n=0:N printf("c%d:%f\n",n,c(n+1)); end; printf("lambda: %f\n",lambda);
f_approx=0; for n=0:N f_approx=f_approx+c(n+1)*s_m.^n; end
printf("normalized sum of square of errors:%f\n",sum((f_approx-f_m).^2)/M);
f_plot=zeros(s_plot); for n=0:N f_plot=f_plot+c(n+1)*s_plot.^n; end
draw_function(s_plot,f_plot,s_m,f_m);
xtitle(title_author,"x","f");
