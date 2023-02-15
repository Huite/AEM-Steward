clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.2.2, Figure 1.20a";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

xmin=-2.5; xmax=2.5; zmin=0; zmax=5.5;
xrange =xmin+(xmax-xmin)*[0:1/500:1];   zrange= zmin+(zmax-zmin)*[0:1/550:1]';   x=(xrange'*ones(zrange'))';    z=(ones(xrange' )*zrange')';
xrangev=xmin+(xmax-xmin)*[1/48:1/24:1]; zrangev=zmin+(zmax-zmin)*[1/56:1/28:1]'; xv=(xrangev'*ones(zrangev'))'; zv=(ones(xrangev')*zrangev')';
pressure_head=zeros(z);
R=.1/365; zlayer=[0,2,3,5,5.5]';
z1=zlayer(1); p1=0;
for ilayer=1:max(size(zlayer))-1
  z0=z1; p0=p1;
  if int(ilayer/2)*2==ilayer then  // for even layers always use fine sand
    Ks=1;    a=5; Fs=Ks/a; k=a/2; 
  else
    Ks=0.02; a=2; Fs=Ks/a; k=a/2;
  end
  c1=(Fs*exp(2*k*p0)-R/(2*k))*exp(2*k*z0);
  c2=R/(2*k);
  z1=zlayer(ilayer+1);       p1=log((c1*exp(-2*k*z1)+c2)/Fs)/(2*k); 
  pressure_head=pressure_head+( log((c1*exp(-2*k*z )+c2)/Fs)/(2*k) ).*(z>z0 & z<=z1);
end
draw_flow_conduction(xrange,zrange,pressure_head,[-2.5:0.1:0]);
xtitle(title_author,"x","z");
draw_segments(xmin*ones(zlayer(2:$-1)),zlayer(2:$-1),xmax*ones(zlayer(2:$-1)),zlayer(2:$-1));
