clear;
clearglobal;

printf("\n  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)\n");
printf("  Oxford University Press.  Supplementary Educational Material.\n\n");
printf("  This work is licensed under the Creative Commons Attribution 4.0 International License.\n");
printf("  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ \n");
printf("  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.\n\n");

title_author="AEM DRSteward Section 1.2.1, Figure 1.17b";
printf("  %s\n\n",title_author);

// execute the file containing functions for visualization
AEM_example_directory=PWD;
exec(AEM_example_directory+"\AEM_DRSteward_functions_visualize.sce");

L1=1000; L2=2000; L4=1000; H3=20; B2=H3; k=30; R=0.0003; n=0.25; hA=28; hD=20;
xmin=0; xmax=L1+L2+L4; zmin=0; zmax=max(hA,hD);
xrange =xmin+(xmax-xmin)*[0:1/501:1];   zrange= zmin+(zmax-zmin)*[0:1/201:1]';   x=(xrange'*ones(zrange'))';    z=(ones(xrange' )*zrange')';
xrangev=xmin+(xmax-xmin)*[1/48:1/24:1]; zrangev=zmin+(zmax-zmin)*[1/56:1/28:1]'; xv=(xrangev'*ones(zrangev'))'; zv=(ones(xrangev')*zrangev')';
hB=sqrt(  R/k*L1*(L2+L4) +(L2+L4)*hA^2/(L1+L2+L4) + (L1)   *hD^2/(L1+L2+L4));
hC=sqrt(  R/k*L4*(L1+L2) +(L4)   *hA^2/(L1+L2+L4) + (L1+L2)*hD^2/(L1+L2+L4));
head=zeros(x);  velocityx=zeros(xv);  velocityz=zeros(xv);

// zone 1
PhiA1=k/2*hA^2; PhiB1=k/2*hB^2;
Phi =-R/2*x .*(x -L1)-(PhiA1-PhiB1)*x /L1+PhiA1;  h =sqrt(2*Phi /k);
Phiv=-R/2*xv.*(xv-L1)-(PhiA1-PhiB1)*xv/L1+PhiA1;  hv=sqrt(2*Phiv/k);
Qx=R*(xv-L1/2)+(PhiA1-PhiB1)/L1;  qx=Qx./hv;                 vx=qx/n;
                                  qz=-(zv./hv).*(R+qx.^2/k); vz=qz/n;
h_zone1 =h.* (0<=x  & x <L1 & z <=h );
vx_zone1=vx.*(0<=xv & xv<L1 & zv<=hv); 
vz_zone1=vz.*(0<=xv & xv<L1 & zv<=hv);

// zone 2
PhiB2=k/2*(hB-B2)^2; PhiC2=k/2*(hC-B2)^2;
Phi =-R/2*(x -L1).*(x -L1-L2)-(PhiB2-PhiC2)*(x -L1)/L2+PhiB2;  h =B2+sqrt(2*Phi /k);  
Phiv=-R/2*(xv-L1).*(xv-L1-L2)-(PhiB2-PhiC2)*(xv-L1)/L2+PhiB2;  hv=B2+sqrt(2*Phiv/k);
Qx=R*(xv-L1-L2/2)+(PhiB2-PhiC2)/L2;  qx=Qx./(hv-B2);                      vx=qx/n;
                                     qz=-((zv-B2)./(hv-B2)).*(R+qx.^2/k); vz=qz/n;
h_zone2 =h.* (L1<=x  & x <L1+L2 & B2<=z  & z <=abs(h) );
vx_zone2=vx.*(L1<=xv & xv<L1+L2 & B2<=zv & zv<=abs(hv));
vz_zone2=vz.*(L1<=xv & xv<L1+L2 & B2<=zv & zv<=abs(hv));

// zone 3
PhiB3=k*H3*hB-k*H3^2/2; PhiC3=k*H3*hC-k*H3^2/2;
Phi=-(PhiB3-PhiC3)*(x-L1)/L2+PhiB3;  h=Phi/(k*H3)+H3/2;
Qx=(PhiB3-PhiC3)/L2*ones(xv); vx=Qx/(H3*n);
                              vz=zeros(vx);
h_zone3 =h.* (L1<=x  & x <L1+L2 & z <H3);
vx_zone3=vx.*(L1<=xv & xv<L1+L2 & zv<H3);
vz_zone3=vz.*(L1<=xv & xv<L1+L2 & zv<H3);

// zone 4
PhiC4=k/2*hC^2; PhiD4=k/2*hD^2;
Phi =-R/2*(x -L1-L2).*(x -L1-L2-L4)-(PhiC4-PhiD4)*(x -L1-L2)/L4+PhiC4;  h =sqrt(2*Phi /k);
Phiv=-R/2*(xv-L1-L2).*(xv-L1-L2-L4)-(PhiC4-PhiD4)*(xv-L1-L2)/L4+PhiC4;  hv=sqrt(2*Phiv/k);
Qx=R*(xv-L1-L2-L4/2)+(PhiC4-PhiD4)/L4;  qx=Qx./hv;                 vx=qx/n;
                                        qz=-(zv./hv).*(R+qx.^2/k); vz=qz/n;;
h_zone4 =h.* (L1+L2<=x  & x <=L1+L2+L4 & z <=h ); 
vx_zone4=vx.*(L1+L2<=xv & xv<=L1+L2+L4 & zv<=hv); 
vz_zone4=vz.*(L1+L2<=xv & xv<=L1+L2+L4 & zv<=hv);

// collect head and velocity from each zone and visualize
head     = h_zone1+ h_zone2+ h_zone3+ h_zone4;
velocityx=vx_zone1+vx_zone2+vx_zone3+vx_zone4;
velocityz=vz_zone1+vz_zone2+vz_zone3+vz_zone4;
draw_flow_conduction(xrange,zrange,real(head),[hD:0.25:hA],[],[],xrangev,zrangev,velocityx,velocityz);
xtitle(title_author,"x","z");
draw_segments([xmin;xmin+L1;xmin+L1+L2;xmin+L1+L2+L4],[zmin;zmin;zmin;zmin],[xmin;xmin+L1;xmin+L1+L2;xmin+L1+L2+L4],[zmax;zmax;zmax;zmax]);
draw_segments([xmin;xmin;xmin+L1],[zmin;zmax;H3],[xmax;xmax;xmin+L1+L2],[zmin;zmax;H3]);
fig=gcf(); fig.figure_size(2)=500; fig.children(2).isoview="off"; fig.visible="on";