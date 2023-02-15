clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 2.2.1, Figure 2.7";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

s_m=[0:1:10]'; s_plot=[0:.1:10]';
f_m=[1.58;1.72;1.70;1.60;0.83;0.46;0.14;0.30;0.51;0.88;1.28]; 
smin=min(s_m); smax=max(s_m); S_m=(s_m-(smin+smax)/2)*2/(smax-smin); S_plot=(s_plot-(smin+smax)/2)*2/(smax-smin);
M=max(size(s_m)); 
N=5;
A=ones(S_m); for n=1:N A=[A,S_m.^n]; end
c=A\f_m;
printf("N:%d\n",N); for n=0:N printf("c%d:%f\n",n,c(n+1)); end
f_approx=0; for n=0:N f_approx=f_approx+c(n+1)*S_m.^n; end
printf("objective:%f\n",sum((f_approx-f_m).^2)/M);
f_plot=zeros(s_plot); for n=0:N f_plot=f_plot+c(n+1)*S_plot.^n; end
draw_function(s_plot,f_plot,s_m,f_m);
xtitle(title_author,"x","f");
