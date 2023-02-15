clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 2.2.2, Figure 2.9";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

s_m=[0:1:10]'; s_plot=[0:.1:10]';
f_m=[1.58;1.72;1.70;1.60;0.83;0.46;0.14;0.30;0.51;0.88;1.28]; 
M=max(size(s_m)); 
theta_m=2*%pi*([1:M]'-1)/M; theta_plot=(s_plot-s_m(1))*2*%pi/(M*(s_m(2)-s_m(1)));
N=3;
A=ones(theta_m); for n=1:N A=[A,cos(n*theta_m),sin(n*theta_m)]; end
c=A\f_m; c0=c(1); ccos=c(2:2:$); csin=c(3:2:$);
printf("N:%d c0:%f\n",N,c0); for n=1:N printf("ccos%d:%f csin%d:%f\n",n,ccos(n),n,csin(n)); end
f_approx=c0*ones(s_m); for n=1:N f_approx=f_approx+ccos(n)*cos(n*theta_m)+csin(n)*sin(n*theta_m); end
printf("objective:%f\n",sum((f_approx-f_m).^2)/M);
f_plot=c0*ones(theta_plot); for n=1:N f_plot=f_plot+ccos(n)*cos(n*theta_plot)+csin(n)*sin(n*theta_plot); end
draw_function(s_plot,f_plot,s_m,f_m);
xtitle(title_author,"x","f");
