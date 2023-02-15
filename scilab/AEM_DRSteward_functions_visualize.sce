//
//  Analytic Element Method: Complex Interactions of Boundaries and Interfaces, David R. Steward (2020)
//  Oxford University Press.  Supplementary Educational Material.
//
//  This work is licensed under the Creative Commons Attribution 4.0 International License.
//  To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/
//  or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.
//
//  AEM_DRSteward_functions_visualize.sce
//

function draw_function(x,f,x_data,f_data)
// x,f: vectors with x-axis locations and values of function to be plotted
// x_data,f_data: vectors with x-axis locations and values of data points to be plotted
  clf(); fig=gcf(); fig.figure_position=[0,0]; fig.figure_size=[720,720];
  drawlater();
  xpoly(x,f);
  ent=gce(); ent.foreground=color("black"); ent.thickness=3; ent.line_style=1; ent.clip_state="on";
  axes=gca();
  axes.data_bounds=[min(x),floor(min(f));max(x),ceil(max(f))]; axes.axes_visible="on"; axes.x_location="bottom";
  axes.thickness=1; axes.grid=[color("forestgreen") color("forestgreen")]; axes.grid_style=[7,7];
  axes.background=color("light yellow"); axes.foreground=color("darkgreen");
  axes.grid_position="foreground";  // move grid in front of plots for newer versions of scilab
  axes.auto_ticks=["on","on","on"]; axes.labels_font_color=color("black");
  axes.font_size=5; axes.font_style=2; axes.sub_tics=[4,4]; axes.title.font_size=5; axes.title.font_style=2; axes.thickness=1;
  deltax=.05*(max(x)-min(x)); deltaf=.05*(ceil(max(f))-floor(min(f)));
  axes.x_label.font_size=5; axes.x_label.font_style=5; axes.x_label.position=[max(x)+deltax,floor(min(f))-deltaf];
  axes.y_label.font_size=5; axes.y_label.font_style=5; axes.y_label.position=[min(x)-deltax,ceil(max(f))+deltaf]; axes.y_label.font_angle=0; 
  drawnow();
  draw_points(x_data,f_data);
endfunction

function draw_flow_conduction(x_function,y_function,phi,phi_interval,psi,psi_interval,x_vector,y_vector,vx,vy)
// x_function: row vector of size (1,Ncol) with x-axis locations where phi and psi have been computed
// y_function: column vector of size (Nrow,1) with y-axis locations where phi and psi have been computed
// phi: real matrix of size(Nrow,Ncol) containing potential function, or [] if this function is not shown
// phi_interval: vector with values of phi to contour, or [] to use default spacing
// psi: real matrix of size(Nrow,Ncol) containing streamfunction, or [] if this function is not shown
// psi_interval: vector with values of psi to contour, or [] to use default spacing
// x_vector: row vector of size (1,Ncol_vector) with x-axis locations where the vector (vx,vy) have been computed
// y_vector: column vector of size (Nrow_vector) with y-axis locations where the vector (vx,vy) have been computed
// vx,vy: matrices of size (Nrow_vector,Ncol_vector) containing vector field component vx and vy, or [] if these functions are not shown
  if phi~=[] then phi=real(phi); end
  if ~exists('phi_interval') then phi_interval=[]; end
  if ~exists('psi')          then psi=[];          end;  if psi~=[] then psi=real(psi); end
  if ~exists('psi_interval') then psi_interval=[]; end
  if ~exists('x_vector')     then x_vector=[];     end
  if ~exists('y_vector')     then y_vector=[];     end
  if ~exists('vx')           then vx=[];           end
  if ~exists('vy')           then vy=[];           end
  if phi==[] then draw_phi=%F; else draw_phi=%T; if phi_interval==[] then phi_interval=[min(phi):(max(phi)-min(phi))/20:max(phi)]; end; end
  if psi==[] then draw_psi=%F; else draw_psi=%T; if psi_interval==[] then psi_interval=[min(psi):(max(psi)-min(psi))/20:max(psi)]; end; end
  if vx==[] | vy==[] then draw_v=%F; else draw_v=%T; end
  
  clf();
  fig=gcf();
  xyratio=(max(x_function)-min(x_function))/(max(y_function)-min(y_function)); 
  ypixel=720; xpixel=ypixel*xyratio; if xpixel>1000 then xpixel=1000; ypixel=xpixel/xyratio; end 
  if draw_phi then xpixel=xpixel/.85; end  // add space for colorbar legend
  fig.figure_position=[0,0]; fig.axes_size=[xpixel,ypixel]; fig.color_map=graycolormap(128);

  if draw_phi then
    drawlater(); 
    if max(size(phi_interval))==1 then n=phi_interval; else n=max(size(phi_interval)); end
    r=0.0+(1.0-0.0)*linspace(0,1,n)';
    g=0.4+(1.0-0.4)*linspace(0,1,n)';
    b=0.0+(0.0-0.0)*linspace(0,1,n)';
    cmap=[r g b];
    fig.color_map=cmap; fig.children.data_bounds=[min(x_function),min(y_function);max(x_function),max(y_function)];
    cmin=min(phi_interval); cmin=round(cmin*10)/10;
    cmax=max(phi_interval); cmax=round(cmax*10)/10;
    colorbar(cmin,cmax);
    fig.children(1).font_size=5; fig.children(1).font_style=2;
    if cmax>=100 then
      axis=fig.children(1); axis.y_ticks.labels=string(int(axis.y_ticks.locations)); axis.axes_bounds(3)=.17; axis.margins(2)=.79;
    end
    contourf(x_function,y_function',phi',phi_interval);
    drawnow();
  end
  
  if draw_psi then
    drawlater(); 
    fig=gcf();
    numentity=max(size(fig.children($).children));
    contour2d(x_function,y_function',psi',psi_interval);
    fig=gcf();
    b2=fig.children($).children(1:max(size(fig.children($).children))-numentity);
    for i=1:1:max(size(b2)) 
      for j=1:1:max(size(b2(i).children))
        if(2*int(j/2)~=j) then
          b2(i).children(j).foreground=color("black");
          b2(i).children(j).thickness=1;
        else
          b2(i).children(j).visible="off";
        end
      end
    end
    drawnow();
  end
  
  if draw_v then
    drawlater(); 
    champ(x_vector,y_vector',vx',vy');
    ent=gce(); ent.arrow_size=1.0;
    drawnow();
  end
 
  if draw_phi then
    drawlater(); 
    fig=gcf();
    numentity=max(size(fig.children($).children));
    xset("fpf"," ");
    contour2d(x_function,y_function',phi',phi_interval);
    fig=gcf();
    b=fig.children($).children(1:max(size(fig.children($).children))-numentity);
    for i=1:1:max(size(b))
      if max(size(b(i).children))==0 then
        c=b(i);
      else
        c=unglue(b(i));
      end
      for j=1:max(size(c))
        c(j).foreground=color("darkgreen"); c(j).thickness=1; c(j).line_style=3;
      end
    end
    drawnow();
  end

  // Adjust grid and labels  
  drawlater(); 
  axes=gca();
  axes.axes_visible="on"; axes.x_location="bottom"; 
  axes.data_bounds=[min(x_function),min(y_function);max(x_function),max(y_function)];
  axes.clip_box=[min(x_function),max(y_function),max(x_function)-min(x_function),max(y_function)-min(y_function)];
  axes.tight_limits="on"; axes.isoview="on"; axes.sub_tics=[4,4]; axes.grid=[color("forestgreen") color("forestgreen")]; axes.grid_style=[7,7];
  axes.background=color("light yellow"); axes.foreground=color("darkgreen");
  axes.grid_position="foreground";  // move grid in front of plots for newer versions of scilab
  axes.auto_ticks=["on","on","on"]; axes.labels_font_color=color("black");
  axes.font_size=5; axes.font_style=2; axes.sub_tics=[4,4]; axes.title.font_size=5; axes.title.font_style=2; axes.thickness=1;
  deltax=.05*(max(x_function)-min(x_function)); deltay=.05*(max(y_function)-min(y_function));
  axes.x_label.font_size=5; axes.x_label.font_style=5; axes.x_label.position=[max(x_function)+deltax,min(y_function)-deltay];
  axes.y_label.font_size=5; axes.y_label.font_style=5; axes.y_label.position=[min(x_function)-deltax,max(y_function)+deltay];
  axes.y_label.font_angle=0; 
  xset("fpf","%.1f");
  drawnow();
endfunction

function draw_waves(x_function,y_function,phi,phi_interval,drawA)
// x_function: row vector of size (1,Ncol) with x-axis locations where phi has been computed
// y_function: column vector of size (Nrow,1) with y-axis locations where phi has been computed
// phi: complex matrix of size(Nrow,Ncol) containing the complex wave function
// phi_interval: vector with values of phi to contour, or number of intervals, or [] to use default spacing
// drawA: %T: draw amplitude of the wave function phi
//        %F: draw cosing of phase for the wave function phi
  if ~exists('phi_interval') then phi_interval=[]; end;  if phi_interval==[] then phi_interval=8; end
  if ~exists('drawA') then drawA=%F; end
  if drawA then
    draw_value=abs(phi);
    if max(size(phi_interval))>1 then
	  draw_interval=phi_interval;
    else
      draw_interval=[0:1/phi_interval:1]*max(abs(phi));
    end
    if max(draw_value)-min(draw_value)<1e-6 then
      printf("\nNot enough change in amplitude to plot\n");
      abort;
    end	
  else
    draw_value=cos(atan(imag(phi),real(phi)));
    if max(size(phi_interval))>1 then
      draw_interval=phi_interval;
    else
      draw_interval=[-1.01,cos(%pi*[-1+1/(2*phi_interval):1/phi_interval:0]),1.01];
    end
  end
  
  clf();
  fig=gcf();
  xyratio=(max(x_function)-min(x_function))/(max(y_function)-min(y_function)); 
  ypixel=720; xpixel=ypixel*xyratio; if xpixel>1000 then xpixel=1000; ypixel=xpixel/xyratio; end 
  xpixel=xpixel/.85; // add space for colorbar legend
  fig.figure_position=[0,0]; fig.axes_size=[xpixel,ypixel]; fig.color_map=graycolormap(128);

  // Draw wave function
  drawlater(); 
  if max(size(draw_interval))==1 then n=draw_interval; else n=max(size(draw_interval)); end
  r=0.0+(1.0-0.0)*linspace(0,1,n)';
  g=0.4+(1.0-0.4)*linspace(0,1,n)';
  b=0.0+(0.0-0.0)*linspace(0,1,n)';
  cmap=[r g b];
  fig.color_map=cmap; fig.children.data_bounds=[min(x_function),min(y_function);max(x_function),max(y_function)];
  cmin=min(draw_interval); cmin=round(cmin*10)/10;
  cmax=max(draw_interval); cmax=round(cmax*10)/10;
  colorbar(cmin,cmax);
  fig.children(1).font_size=5; fig.children(1).font_style=2;
  if cmax>=100 then
    axis=fig.children(1); axis.y_ticks.labels=string(int(axis.y_ticks.locations)); axis.axes_bounds(3)=.17; axis.margins(2)=.79;
  end
  contourf(x_function,y_function',draw_value',draw_interval);
  drawnow();
 
  // Draw contours
  drawlater(); 
  fig=gcf();
  numentity=max(size(fig.children($).children));
  xset("fpf"," ");
  contour2d(x_function,y_function',draw_value',draw_interval);
  fig=gcf();
  b=fig.children($).children(1:max(size(fig.children($).children))-numentity);
  for i=1:1:max(size(b))
    if max(size(b(i).children))==0 then
      c=b(i);
    else
      c=unglue(b(i));
    end
    for j=1:max(size(c))
      c(j).foreground=color("darkgreen"); c(j).thickness=1; c(j).line_style=1;
    end
  end
  drawnow();
 
  if drawA then // draw thicker line for A=1;
    drawlater();
    a=gcf();
    numentity=max(size(a.children($).children));
    xset("fpf"," ");
    contour2d(x_function,y_function',draw_value',[0 1]);
    a=gcf();
    b=a.children($).children(1:max(size(a.children($).children))-numentity);
    for i=1:1:max(size(b))
      if max(size(b(i).children))==0 then
        c=b(i);
      else
        c=unglue(b(i));
      end
      for j=1:max(size(c))
        c(j).foreground=color("black");
        c(j).thickness=3;
      end
    end
    drawnow();
  end

  // Adjust grid and labels  
  drawlater(); 
  axes=gca();
  axes.axes_visible="on"; axes.x_location="bottom"; 
  axes.data_bounds=[min(x_function),min(y_function);max(x_function),max(y_function)];
  axes.clip_box=[min(x_function),max(y_function),max(x_function)-min(x_function),max(y_function)-min(y_function)];
  axes.tight_limits="on"; axes.isoview="on"; axes.sub_tics=[4,4]; axes.grid=[color("forestgreen") color("forestgreen")]; axes.grid_style=[7,7];
  axes.background=color("light yellow"); axes.foreground=color("darkgreen");
  axes.grid_position="foreground";  // move grid in front of plots for newer versions of scilab
  axes.auto_ticks=["on","on","on"]; axes.labels_font_color=color("black");
  axes.font_size=5; axes.font_style=2; axes.sub_tics=[4,4]; axes.title.font_size=5; axes.title.font_style=2; axes.thickness=1;
  deltax=.05*(max(x_function)-min(x_function)); deltay=.05*(max(y_function)-min(y_function));
  axes.x_label.font_size=5; axes.x_label.font_style=5; axes.x_label.position=[max(x_function)+deltax,min(y_function)-deltay];
  axes.y_label.font_size=5; axes.y_label.font_style=5; axes.y_label.position=[min(x_function)-deltax,max(y_function)+deltay];
  axes.y_label.font_angle=0; 
  xset("fpf","%.1f");
  drawnow();
endfunction

function draw_stress(x_function,y_function,sigmax,sigmay,tauxy,stress_interval,x_vector,y_vector,vec,drawstress)
// x_function: row vector of size (1,Ncol) with x-axis locations where stress has been computed
// y_function: column vector of size (Nrow,1) with y-axis locations where stress has been computed
// sigmax: real matrix of size(Nrow,Ncol) containing this stress component
// sigmay: real matrix of size(Nrow,Ncol) containing this stress component
// tauxy: real matrix of size(Nrow,Ncol) containing this stress component
// stress_interval: vector with values of mean stress or maximum shear to contour, or [] to use default spacing
// x_vector: row vector of size (1,Ncol_vector) with x-axis locations where the stress or displacement have been computed
// y_vector: column vector of size (Nrow_vector) with y-axis locations where the stress or displacement have been computed
// vec: vector of matrices of size (Nrow_vector,Ncol_vector) containing stress components or displacement components corresponding to drawstress
// drawstress %T: plot mean stress (isostress) and show principal stresses using vec.sigmax,vec.sigmay,vec.tauxy at x_vector,y_vector locations
//            %F: plot max shear   (isoshear)  and show displacement       using vec.ux,vec.uy                   at x_vector,y_vector locations
  if ~exists('drawstress') then drawstress=%T; end
  
  clf();
  fig=gcf();
  xyratio=(max(x_function)-min(x_function))/(max(y_function)-min(y_function)); 
  ypixel=720; xpixel=ypixel*xyratio; if xpixel>1000 then xpixel=1000; ypixel=xpixel/xyratio; end 
  xpixel=xpixel/.85; // add space for colorbar legend
  fig.figure_position=[0,0]; fig.axes_size=[xpixel,ypixel]; fig.color_map=graycolormap(128);

  if drawstress then
    // isostress lines of uniform mean stress and principal stresses
    sigma1=(sigmax+sigmay)/2+(((sigmax-sigmay)/2).^2+tauxy.^2).^(0.5);
    sigma2=(sigmax+sigmay)/2-(((sigmax-sigmay)/2).^2+tauxy.^2).^(0.5);
    sigma_isostress=(sigma1+sigma2)/2;  // iscopachs are lines of of constant mean stress
    if stress_interval==[] then
      sigma_isostress_min=min(sigma_isostress); sigma_isostress_max=max(sigma_isostress); 
      sigma_isostress_range=0.5*(sigma_isostress_min+sigma_isostress_max)+1.5*[-1:.05:1]*(0.5*(sigma_isostress_max-sigma_isostress_min));;
    else
      sigma_isostress_range=stress_interval;
    end
  
     drawlater(); 
    if max(size(sigma_isostress_range))==1 then n=sigma_isostress_range; else n=max(size(sigma_isostress_range)); end
    r=0.0+(1.0-0.0)*linspace(0,1,n)';
    g=0.4+(1.0-0.4)*linspace(0,1,n)';
    b=0.0+(0.0-0.0)*linspace(0,1,n)';
    cmap=[r g b];
    fig.color_map=cmap; fig.children.data_bounds=[min(x_function),min(y_function);max(x_function),max(y_function)];
    cmin=min(sigma_isostress_range); cmin=round(cmin*10)/10;
    cmax=max(sigma_isostress_range); cmax=round(cmax*10)/10;
    colorbar(cmin,cmax);
    fig.children(1).font_size=5; fig.children(1).font_style=2;
    if cmax>=100 then
      axis=fig.children(1); axis.y_ticks.labels=string(int(axis.y_ticks.locations)); axis.axes_bounds(3)=.17; axis.margins(2)=.79;
    end
    contourf(x_function,y_function',sigma_isostress',sigma_isostress_range);
    drawnow();

    drawlater(); 
    fig=gcf();
    numentity=max(size(fig.children($).children));
    xset("fpf"," ");
    contour2d(x_function,y_function',sigma_isostress',sigma_isostress_range);
    fig=gcf();
    b=fig.children($).children(1:max(size(fig.children($).children))-numentity);
    for i=1:1:max(size(b))
      if max(size(b(i).children))==0 then
        c=b(i);
      else
        c=unglue(b(i));
      end
      for j=1:max(size(c))
        c(j).foreground=color("darkgreen"); c(j).thickness=1; c(j).line_style=1;
      end
    end
    drawnow();

    // compute eigenvalues (principal stresses) and eigenvectors (principal directions) for vec.sigmax, vec.sigmay, vec.tauxy
    vec.sigma1=(vec.sigmax+vec.sigmay)/2+(((vec.sigmax-vec.sigmay)/2).^2+vec.tauxy.^2).^(0.5);
    vec.sigma2=(vec.sigmax+vec.sigmay)/2-(((vec.sigmax-vec.sigmay)/2).^2+vec.tauxy.^2).^(0.5);
    vec.taumax=                          (((vec.sigmax-vec.sigmay)/2).^2+vec.tauxy.^2).^(0.5);
    vec.eigenvalue1=zeros(vec.sigmax);  vec.eigenvector1_dx=zeros(vec.sigmax);  vec.eigenvector1_dy=zeros(vec.sigmax); 
    vec.eigenvalue2=zeros(vec.sigmax);  vec.eigenvector2_dx=zeros(vec.sigmax);  vec.eigenvector2_dy=zeros(vec.sigmax); 
    for irow=1:size(vec.sigmax,1)
      for jcol=1:size(vec.sigmax,2)
         [eigenvector,eigenvalue]=spec([vec.sigmax(irow,jcol),vec.tauxy(irow,jcol);vec.tauxy(irow,jcol),vec.sigmay(irow,jcol)])
         vec.eigenvalue2(irow,jcol)=eigenvalue(1,1); vec.eigenvector2_dx(irow,jcol)=eigenvector(1,1); vec.eigenvector2_dy(irow,jcol)=eigenvector(2,1)
         vec.eigenvalue1(irow,jcol)=eigenvalue(2,2); vec.eigenvector1_dx(irow,jcol)=eigenvector(1,2); vec.eigenvector1_dy(irow,jcol)=eigenvector(2,2)
      end
    end
    // compute straight lines in principal directions scaled by principal stress and of length that they will not overlap
    deltax=0.5*abs(x_vector(2)-x_vector(1))  // half distance between grid points x-direction
    deltay=0.5*abs(y_vector(2)-y_vector(1))  // half distance between grid points y-direction
    dsigmax=max( max(abs(vec.eigenvalue1.*vec.eigenvector1_dx)) , max(abs(vec.eigenvalue2.*vec.eigenvector2_dx)));
    if dsigmax>1e-8 then
      scalex=deltax/dsigmax;
    else
      scalex=deltax;
    end
    dsigmay=max( max(abs(vec.eigenvalue1.*vec.eigenvector1_dy)) , max(abs(vec.eigenvalue2.*vec.eigenvector2_dy)));
    if dsigmay>1e-8 then
      scaley=deltay/dsigmay;
    else
      scaley=deltay;
    end
    scale=min(scalex,scaley); scale=scale*.80;

	// draw principal stresses
    drawlater();
    sigmacolor=color(42,13,70)+(color("white")-color(42,13,70))*(vec.eigenvalue1>0);
    sigmacolor=matrix(sigmacolor,1,-1);
    vxstart=(x_vector'*ones(y_vector'))';  vystart=(ones(x_vector')*y_vector')';
    vxend=vxstart+scale*vec.eigenvalue1.*vec.eigenvector1_dx;
    vyend=vystart+scale*vec.eigenvalue1.*vec.eigenvector1_dy;
    vx=[matrix(vxstart,1,-1);matrix(vxend,1,-1)]; vy=[matrix(vystart,1,-1);matrix(vyend,1,-1)]; 
    xarrows(vx,vy);  e=gce(); e.arrow_size=0; e.thickness=2; e.segs_color=ones(e.segs_color).*sigmacolor;
    vxend=vxstart-scale*vec.eigenvalue1.*vec.eigenvector1_dx;
    vyend=vystart-scale*vec.eigenvalue1.*vec.eigenvector1_dy;
    vx=[matrix(vxstart,1,-1);matrix(vxend,1,-1)]; vy=[matrix(vystart,1,-1);matrix(vyend,1,-1)]; 
    xarrows(vx,vy);  e=gce(); e.arrow_size=0; e.thickness=2; e.segs_color=ones(e.segs_color).*sigmacolor;
    sigmacolor=color(42,13,70)+(color("white")-color(42,13,70))*(vec.eigenvalue2>0);
    sigmacolor=matrix(sigmacolor,1,-1);
    vxend=vxstart+scale*vec.eigenvalue2.*vec.eigenvector2_dx;
    vyend=vystart+scale*vec.eigenvalue2.*vec.eigenvector2_dy;
    vx=[matrix(vxstart,1,-1);matrix(vxend,1,-1)]; vy=[matrix(vystart,1,-1);matrix(vyend,1,-1)]; 
    xarrows(vx,vy);  e=gce(); e.arrow_size=0; e.thickness=2; e.segs_color=ones(e.segs_color).*sigmacolor;
    vxend=vxstart-scale*vec.eigenvalue2.*vec.eigenvector2_dx;
    vyend=vystart-scale*vec.eigenvalue2.*vec.eigenvector2_dy;
    vx=[matrix(vxstart,1,-1);matrix(vxend,1,-1)]; vy=[matrix(vystart,1,-1);matrix(vyend,1,-1)]; 
    xarrows(vx,vy);  e=gce(); e.arrow_size=0; e.thickness=2; e.segs_color=ones(e.segs_color).*sigmacolor;
    drawnow();
  else
    // isochromatic lines of uniform shear stress (isoshear) and displacement vector
    taumax=(((sigmax-sigmay)/2).^2+tauxy.^2).^(0.5);
    if stress_interval==[] then
      taumax_isoshear_min=min(taumax); taumax_isoshear_max=max(taumax); 
      taumax_isoshear_range=0.5*(taumax_isoshear_min+taumax_isoshear_max)+1.5*[-1:.1:1]*(0.5*(taumax_isoshear_max-taumax_isoshear_min));
    else
      taumax_isoshear_range=stress_interval;
    end
 
    drawlater(); 
    if max(size(taumax_isoshear_range))==1 then n=taumax_isoshear_range; else n=max(size(taumax_isoshear_range)); end
    r=0.0+(1.0-0.0)*linspace(0,1,n)';
    g=0.4+(1.0-0.4)*linspace(0,1,n)';
    b=0.0+(0.0-0.0)*linspace(0,1,n)';
    cmap=[r g b];
    fig.color_map=cmap; fig.children.data_bounds=[min(x_function),min(y_function);max(x_function),max(y_function)];
    cmin=min(taumax_isoshear_range); cmin=round(cmin*10)/10;
    cmax=max(taumax_isoshear_range); cmax=round(cmax*10)/10;
    colorbar(cmin,cmax);
    fig.children(1).font_size=5; fig.children(1).font_style=2;
    if cmax>=100 then
      axis=fig.children(1); axis.y_ticks.labels=string(int(axis.y_ticks.locations)); axis.axes_bounds(3)=.17; axis.margins(2)=.79;
    end
    contourf(x_function,y_function',taumax',taumax_isoshear_range);
    drawnow();

    drawlater(); 
    fig=gcf();
    numentity=max(size(fig.children($).children));
    xset("fpf"," ");
    contour2d(x_function,y_function',taumax',taumax_isoshear_range);
    fig=gcf();
    b=fig.children($).children(1:max(size(fig.children($).children))-numentity);
    for i=1:1:max(size(b))
      if max(size(b(i).children))==0 then
        c=b(i);
      else
        c=unglue(b(i));
      end
      for j=1:max(size(c))
        c(j).foreground=color("darkgreen"); c(j).thickness=1; c(j).line_style=2;
      end
    end
    drawnow();

	// visualize displacement
    drawlater(); 
    champ(x_vector,y_vector',vec.ux',vec.uy');
    ent=gce(); ent.arrow_size=1.0;
    drawnow();
  end

  // Adjust grid and labels  
  drawlater(); 
  axes=gca();
  axes.axes_visible="on"; axes.x_location="bottom"; 
  axes.data_bounds=[min(x_function),min(y_function);max(x_function),max(y_function)];
  axes.clip_box=[min(x_function),max(y_function),max(x_function)-min(x_function),max(y_function)-min(y_function)];
  axes.tight_limits="on"; axes.isoview="on"; axes.sub_tics=[4,4]; axes.grid=[color("forestgreen") color("forestgreen")]; axes.grid_style=[7,7];
  axes.background=color("light yellow"); axes.foreground=color("darkgreen");
  axes.grid_position="foreground";  // move grid in front of plots for newer versions of scilab
  axes.auto_ticks=["on","on","on"]; axes.labels_font_color=color("black");
  axes.font_size=5; axes.font_style=2; axes.title.font_size=5; axes.title.font_style=2; axes.thickness=1;
  deltax=.05*(max(x_function)-min(x_function)); deltay=.05*(max(y_function)-min(y_function));
  axes.x_label.font_size=5; axes.x_label.font_style=5; axes.x_label.position=[max(x_function)+deltax,min(y_function)-deltay];
  axes.y_label.font_size=5; axes.y_label.font_style=5; axes.y_label.position=[min(x_function)-deltax,max(y_function)+deltay];
  axes.y_label.font_angle=0; 
  xset("fpf","%.1f");
  drawnow();
endfunction

function draw_segments(xend1,yend1,xend2,yend2)
// xend1,yend1: vectors containing locations of first endpoint of each line segment
// xend2,yend2: vectors containing locations of second endpoint of each line segment
  drawlater();
  for i=1:max(size(xend1))
    plot([xend1(i),xend2(i)],[yend1(i),yend2(i)]);
    ent=gce(); ent.children.thickness=3; ent.children.foreground=color(42,13,70);
  end
  drawnow();
endfunction

function draw_points(x,y)
// x,y: vectors containing locations of points
  drawlater();
  for i=1:max(size(x))
    xpoly(x(i),y(i));
    ent=gce(); ent.line_mode="off"; ent.mark_mode="on"; ent.mark_size_unit="point"; ent.mark_foreground=color(42,13,70); ent.clip_state="clipgrf";
    ent.mark_size=7; ent.thickness=2; ent.mark_style=10; // (0: filled circle, 1:plus, 4:filled diamond, 10:star)
  end
  drawnow();
endfunction
