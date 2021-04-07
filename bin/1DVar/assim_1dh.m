function [posterior,diagnostics,representers]=assim_1dh(prior,obs,verb,bkgd)
%
% [posterior,diagnostics]=assim_1dh(prior,obs)
%
% Assimilate data into 1DH wave and current model.
%
% INPUTS:
%
% prior.x = model grid
% prior.h = bathymetry
% prior.theta0, prior.H0 = boundary conditions
% prior.ka_drag = bed roughness
% prior.tauw = wind stress in m2/s2 units (default = 0 if not included)
%
% prior.Ch, prior.Ctheta0, prior.CH0, prior.Cka = covariance matrices for h,
% theta0, H0, ka_drag
%
% obs.H.ind = model indeces for rms wave height observations
% obs.H.d = rms wave height data
% obs.H.e = error stdev for rms wave height
% obs.v.* = as in obs.H, but for longshore current
%
% NOTE, if any of obs.* are missing or empty, they will be ignored
%
% OUTPUTS:
%
% posterior = struct with same fields as prior but updated by assimilation
%
% diagnostics = struct with diagnostic fields for analyzing the update.  See
% comments near end of this code for descriptions
%
addpath '/home/wilsongr/work/funded_projects/USCRP_DA/project/cmtb/bin/1DVar/waveModel/';
warning('off','all');
if(~exist('verb'))
  verb=1;
end

hmin=.2;
nitermax=500;

% unpack some variables
Ch=prior.Ch;
CH0=prior.CH0;
Ctheta0=prior.Ctheta0;
Cka=prior.Cka;
x=prior.x;
nx=length(x);
nt=length(prior);

% if any obs types missing, set to empty
fld='Hv';
noobs.ind=[];
noobs.d=[];
noobs.e=[];
for i=1:length(fld)
  if(~isfield(obs,fld(i)))
    obs=setfield(obs,fld(i),noobs);
  end
end

% Adjust prior tauw to account for possible pressure gradients.  The eqn for
% v is copied from waveModel.m and directly inverted to solve for the force
% Fy, using alongshore current observations on the shelf (x>800m).  This Fy
% is then assigned to tauw=Fy, replacing the input wind stress.  It should
% be noted this total force Fy includes other contributions other than wind
% stress (dp/dy); so assigning it to tauw is an abuse of notation, but
% avoids the need for other coding changes such as adding a separate dp/dy
% term to waveModel.m.
imin=min(find(x>800));
ind=find(obs.v.ind>imin);
if(isempty(ind))
  prior.tauw = obs.tauw;  % default: if vshelf observations are unavailable, just use wind stress
else
  vshelf = obs.v.d(ind);
  h = interp1(x,prior.h,x(obs.v.ind(ind)));
  g=9.8126;
  for i=1:length(ind)
    k(i)=fzero(@(k)prior.sigma^2-g*k.*tanh(k.*h(i)),prior.sigma./sqrt(g*h(i)),optimset('Display','off'));
  end
  a=1.16;  % empirical constant
  Cd=0.015*(prior.ka_drag./h).^(1/3);
  urms=1.416*prior.H0.*prior.sigma./(4*sinh(k.*h));
  % v2 = sqrt( (a*Cd.*urms).^4 + 4*(Cd.*Fy).^2 )./(2*Cd.^2) - (a*urms).^2/2;  % from waveModel.m
  Fy = sqrt( ( ( ( vshelf.^2 + (a*urms).^2/2 ).*(2*Cd.^2) ).^2 - (a*Cd.*urms).^4 )./(4*Cd.^2) );
  Fy = -abs(Fy).*sign(mean(vshelf));
  prior.tauw = mean(Fy);  % override user-input tauw with inverted version
  clear h k a Cd urms
  % ind=setdiff(1:length(obs.v.ind),ind);  % optional: toss out obs used for calculating Fy above
  % for fld={'ind','d','e'}
  %   fld=cell2mat(fld);
  %   this=getfield(obs.v,fld);
  %   obs.v=setfield(obs.v,fld,this(ind));
  % end
end

% initialize prior model state using NL model
if(verb)
  disp('running prior forward model...')
end
prior=waveModel(x,prior.H0,prior.theta0,prior);
if(verb)
  disp('...done')
end

% outer loop
eps=nan;
for n=1:nitermax
  disp(['iteration ' num2str(n) ', itermax = ' num2str(nitermax) ', eps = ' num2str(eps)])
  % keyboard

  % define background and forecast states for this outer loop iteration
  if(n==1)
    if(~exist('bkgd'))
      bkgd=prior;
    else
      disp('using provided bkgd state')
    end
    fcst=prior;
  else
    bkgd=posterior;
    fcst=prior;
    tl_h=fcst.h-bkgd.h;
    tl_H0=fcst.H0-bkgd.H0;
    tl_theta0=fcst.theta0-bkgd.theta0;
    tl_ka_drag=fcst.ka_drag-bkgd.ka_drag;
    [tl_H,tl_theta,tl_v]=tl_waveModel(x,tl_h,tl_H0,tl_theta0,tl_ka_drag,bkgd);
    fcst.H=bkgd.H+tl_H;
    fcst.theta=bkgd.theta+tl_theta;
    fcst.v=bkgd.v+tl_v;
  end

  % initialize representer sub-matrices
  %
  % LEGEND:
  %
  % R_XY = sensitivity of observation type Y, to delta-perturbations of
  % observation X.  These are sub-blocks of L*M*C*M'*L' (note, C is the
  % prior covariance).
  %
  % r_X = sensitivity of model input vector phi = [h; H0; theta0; ka_drag])
  % to delta-perturbations of observation X.  These are rows of C*M'*L'.
  %
  % r_XY = sensitivity of model variable Y to delta-perturbations of
  % observation X.  This is used for diagnostic purposes only, not for
  % assimilation, and is output in struct 'rdiag'.  These are just like r_X
  % (above) but are for variables that aren't treated as model inputs
  % (v,H,k,etc.).
  %
  % ad_X_Y = sensitivity of model input variable Y (one of the phi
  % variables) to delta-perturbations of observation X.  These are direct
  % outputs of the adjoint model, i.e. M'*L', with no smoothing by the prior
  % covariance.  Stored for diagnostic purposes only.
  %
  clear R_* r_* rep_*
  r_H=zeros(nx+3,length(obs.H.ind));
  r_v=zeros(nx+3,length(obs.v.ind));
  ad_H_h      =zeros(length(obs.H.ind),nx);
  ad_H_H0     =zeros(length(obs.H.ind),1);
  ad_H_theta0 =zeros(length(obs.H.ind),1);
  ad_H_ka_drag=zeros(length(obs.H.ind),1);
  ad_v_h      =zeros(length(obs.v.ind),nx);
  ad_v_H0     =zeros(length(obs.v.ind),1);
  ad_v_theta0 =zeros(length(obs.v.ind),1);
  ad_v_ka_drag=zeros(length(obs.v.ind),1);
  R_HH=zeros(length(obs.H.ind),length(obs.H.ind));
  R_Hv=zeros(length(obs.H.ind),length(obs.v.ind));
  R_vH=zeros(length(obs.v.ind),length(obs.H.ind));
  R_vv=zeros(length(obs.v.ind),length(obs.v.ind));
  r_HH=zeros(length(obs.H.ind),nx);
  r_Hv=zeros(length(obs.H.ind),nx);
  r_vH=zeros(length(obs.v.ind),nx);
  r_vv=zeros(length(obs.v.ind),nx);

  % compute representers for wave height: apply delta-perturbations to
  % observations of type X, to compute (a) sensitivity of model inputs psi
  % (r_* matrices), and (b) sensitivity of observations of type Y (R_*
  % matrices)
  for i=1:length(obs.H.ind)
    comb=zeros(nx,1);
    comb(obs.H.ind(i))=1;  % data functional (delta-fn, aka identity matrix)
    [ad_h,ad_H0,ad_theta0,ad_ka_drag]=ad_waveModel(x,comb,0*comb,0*comb,0*comb,bkgd);  % I*M', where M is the TL model and I is identity
    r_H(:,i)=[Ch*ad_h; CH0*ad_H0; Ctheta0*ad_theta0; Cka*ad_ka_drag];  % = C*M'
    ad_H_h(i,:)=ad_h;
    ad_H_H0(i,:)=ad_H0;
    ad_H_theta0(i,:)=ad_theta0;
    ad_H_ka_drag(i,:)=ad_ka_drag;
    [tl_H,tl_a,tl_v,tl_k]=tl_waveModel(x,Ch*ad_h,CH0*ad_H0,Ctheta0*ad_theta0,Cka*ad_ka_drag,bkgd);
    r_HH(i,:)=tl_H;
    r_Hv(i,:)=tl_v;
    r_Hh(i,:)=Ch*ad_h;
    R_HH(i,:)=tl_H(obs.H.ind);
    R_Hv(i,:)=tl_v(obs.v.ind);
  end

  % repeat to compute representers for longshore current
  for i=1:length(obs.v.ind)
    comb=zeros(nx,1);
    comb(obs.v.ind(i))=1;  % data functional (delta-fn, aka identity matrix)
    [ad_h,ad_H0,ad_theta0,ad_ka_drag]=ad_waveModel(x,0*comb,0*comb,comb,0*comb,bkgd);  % I*M', where M is the TL model and I is identity
    r_v(:,i)=[Ch*ad_h; CH0*ad_H0; Ctheta0*ad_theta0; Cka*ad_ka_drag];  % = C*M'
    ad_v_h(i,:)=ad_h;
    ad_v_H0(i,:)=ad_H0;
    ad_v_theta0(i,:)=ad_theta0;
    ad_v_ka_drag(i,:)=ad_ka_drag;
    [tl_H,tl_a,tl_v,tl_k]=tl_waveModel(x,Ch*ad_h,CH0*ad_H0,Ctheta0*ad_theta0,Cka*ad_ka_drag,bkgd);
    r_vH(i,:)=tl_H;
    r_vv(i,:)=tl_v;
    r_vh(i,:)=Ch*ad_h;
    R_vH(i,:)=tl_H(obs.H.ind);
    R_vv(i,:)=tl_v(obs.v.ind);
  end

  % assemble matrices for updating
  Cd=diag([obs.H.e(:);
           obs.v.e(:)].^2);
  CMt=[r_H r_v];
  d=[obs.H.d(:);
     obs.v.d(:)];
  Lu=[fcst.H(obs.H.ind);
      fcst.v(obs.v.ind)];
  R=[R_HH R_Hv;
     R_vH R_vv];
  obstype=[repmat('H',[length(obs.H.d) 1]);
           repmat('v',[length(obs.v.d) 1])];

  % compute the update
  hedge=1; %tanh(n/5);  % reduce magnitude of update for 1st few iterations
  update=hedge*CMt*inv(R+Cd)*(d-Lu);
  if(n>1)
    update=.5*(u0+update);
  end
  u0=update;
  posterior=prior;
  posterior.H0=prior.H0+update(nx+1);
  posterior.theta0=prior.theta0+update(nx+2);
  posterior.ka_drag=max(1e-4,prior.ka_drag+update(nx+3));
  posterior.h = prior.h + update(1:nx);
  posterior.h(posterior.h<hmin)=hmin;
  C2=blkdiag(Ch,CH0,Ctheta0,Cka)-CMt*inv(R+Cd)*CMt';
  v2=diag(C2);
  posterior.hErr=sqrt(max(0,v2(1:nx)));
  posterior.Ch=C2(1:nx,1:nx);
  posterior.CH0=v2(nx+2);
  posterior.Ctheta0=v2(nx+2);
  posterior.Cka=v2(nx+3);
  posterior=waveModel(x,posterior.H0,posterior.theta0,posterior);

  % show the results
  if(verb)
    %disp('plot')
    clf
    subplot(221), hold on
    plot(prior.x,prior.h,'r')
    plot(prior.x,prior.h-sqrt(diag(prior.Ch)),'r--')
    plot(prior.x,prior.h+sqrt(diag(prior.Ch)),'r--')
    plot(posterior.x,posterior.h,'b')
    plot(posterior.x,posterior.h-posterior.hErr,'b--')
    plot(posterior.x,posterior.h+posterior.hErr,'b--')
    set(gca,'ydir','reverse')
    ylabel('h [m]')
    % legend('prior','inverted')
    subplot(222), hold on
    plot(prior.x,prior.H,'r')
    plot(posterior.x,posterior.H,'b')
    errorbar(x(obs.H.ind),obs.H.d,obs.H.e,'ko')
    ylabel('H_{rms} [m]')
    title(['H0 prior=' num2str(prior.H0,2) 'm, final=' num2str(posterior.H0,2) 'm']);
    subplot(223), hold on
    plot(prior.x,rad2deg(prior.theta),'r')
    plot(posterior.x,rad2deg(posterior.theta),'b')
    ylabel('theta [deg]')
    title(['theta0 prior=' num2str(rad2deg(prior.theta0),2) 'deg, final=' num2str(rad2deg(posterior.theta0),2) 'deg']);
    subplot(224), hold on
    plot(prior.x,prior.v,'r')
    plot(posterior.x,posterior.v,'b')
    errorbar(x(obs.v.ind),obs.v.d,obs.v.e,'ko')
    ylabel('v [m/s]')
    title(['ka prior=' num2str(prior.ka_drag,3) 'm, final=' num2str(posterior.ka_drag,3) 'm']);
    for i=1:4
      subplot(2,2,i)
      box on
    end
    %disp('plot2')
    %pause(.1)
  end

  % check for convergence
  eps=sum((posterior.h-bkgd.h).^2)/trace(Ch);
  %if(CH0>0)
  %  eps=eps+(posterior.H0-bkgd.H0)^2/CH0;
  %end
  %if(Ctheta0>0)
  %  eps=eps+(posterior.theta0-bkgd.theta0)^2/Ctheta0;
  %end
  %if(Cka>0)
  %  eps=eps+(posterior.ka_drag-bkgd.ka_drag)^2/Cka;
  %end
  if(eps<1e-4)
      disp(eps)
    break;
  end

end  % outer loop iterations

% OPTIONAL: for diagnostics, repeat the representer calculations for both
% variables, but this time do ALL the gridpoints.  These can be used for
% prior and posterior covariances for v, H.  This is not used for the outer
% loop updates, so only need to do it once at the end
if(nargout>2)
  disp('Calculating full representer matrix as requested (nargout>2)')

  adj_H_h      =zeros(nx,nx);
  adj_H_H0     =zeros(nx,1);
  adj_H_theta0 =zeros(nx,1);
  adj_H_ka_drag=zeros(nx,1);
  adj_v_h      =zeros(nx,nx);
  adj_v_H0     =zeros(nx,1);
  adj_v_theta0 =zeros(nx,1);
  adj_v_ka_drag=zeros(nx,1);
  rep_HH=zeros(nx);
  rep_Hv=zeros(nx);
  rep_Hh=zeros(nx);
  rep_vH=zeros(nx);
  rep_vv=zeros(nx);
  rep_vh=zeros(nx);
  rep_HH_post=zeros(nx);
  rep_Hv_post=zeros(nx);
  rep_Hh_post=zeros(nx);
  rep_vH_post=zeros(nx);
  rep_vv_post=zeros(nx);
  rep_vh_post=zeros(nx);

  % wave height representers
  for i=1:nx
    comb=zeros(nx,1);
    comb(i)=1;
    [ad_h,ad_H0,ad_theta0,ad_ka_drag]=ad_waveModel(x,comb,0*comb,0*comb,0*comb,bkgd);  % I*M', where M is the TL model and I is identity
    adj_H_h(i,:)=ad_h;
    adj_H_H0(i,:)=ad_H0;
    adj_H_theta0(i,:)=ad_theta0;
    adj_H_ka_drag(i,:)=ad_ka_drag;
    [tl_H,tl_a,tl_v,tl_k]=tl_waveModel(x,Ch*ad_h,CH0*ad_H0,Ctheta0*ad_theta0,Cka*ad_ka_drag,bkgd);
    rep_HH(i,:)=tl_H;
    rep_Hv(i,:)=tl_v;
    rep_Hh(i,:)=Ch*ad_h;
    [tl_H,tl_a,tl_v,tl_k]=tl_waveModel(x,posterior.Ch*ad_h,posterior.CH0*ad_H0,posterior.Ctheta0*ad_theta0,posterior.Cka*ad_ka_drag,bkgd);
    rep_HH_post(i,:)=tl_H;
    rep_Hv_post(i,:)=tl_v;
    rep_Hh_post(i,:)=posterior.Ch*ad_h;
  end

  % longshore current representers
  for i=1:nx
    comb=zeros(nx,1);
    comb(i)=1;  % data functional (delta-fn, aka identity matrix)
    [ad_h,ad_H0,ad_theta0,ad_ka_drag]=ad_waveModel(x,0*comb,0*comb,comb,0*comb,bkgd);  % I*M', where M is the TL model and I is identity
    adj_v_h(i,:)=ad_h;
    adj_v_H0(i,:)=ad_H0;
    adj_v_theta0(i,:)=ad_theta0;
    adj_v_ka_drag(i,:)=ad_ka_drag;
    [tl_H,tl_a,tl_v,tl_k]=tl_waveModel(x,Ch*ad_h,CH0*ad_H0,Ctheta0*ad_theta0,Cka*ad_ka_drag,bkgd);
    rep_vH(i,:)=tl_H;
    rep_vv(i,:)=tl_v;
    rep_vh(i,:)=Ch*ad_h;
    [tl_H,tl_a,tl_v,tl_k]=tl_waveModel(x,posterior.Ch*ad_h,posterior.CH0*ad_H0,posterior.Ctheta0*ad_theta0,posterior.Cka*ad_ka_drag,bkgd);
    rep_vH_post(i,:)=tl_H;
    rep_vv_post(i,:)=tl_v;
    rep_vh_post(i,:)=posterior.Ch*ad_h;
  end

  representers.adj_H_h      =adj_H_h      ;
  representers.adj_H_H0     =adj_H_H0     ;
  representers.adj_H_theta0 =adj_H_theta0 ;
  representers.adj_H_ka_drag=adj_H_ka_drag;
  representers.adj_v_h      =adj_v_h      ;
  representers.adj_v_H0     =adj_v_H0     ;
  representers.adj_v_theta0 =adj_v_theta0 ;
  representers.adj_v_ka_drag=adj_v_ka_drag;
  representers.H_H=rep_HH;
  representers.H_v=rep_Hv;
  representers.H_h=rep_Hh;
  representers.v_H=rep_vH;
  representers.v_v=rep_vv;
  representers.v_h=rep_vh;
  representers.H_H_post=rep_HH_post;
  representers.H_v_post=rep_Hv_post;
  representers.H_h_post=rep_Hh_post;
  representers.v_H_post=rep_vH_post;
  representers.v_v_post=rep_vv_post;
  representers.v_h_post=rep_vh_post;

end

%---------------------------------------
% remaining code: reformat output variables in diagnostics struct, for
% convenience
%---------------------------------------
clear diagnostics

diagnostics.niter=n;
diagnostics.eps=eps;

% terms in update equation
diagnostics.update.CMt=CMt;
diagnostics.update.R  =R  ;
diagnostics.update.Cd =Cd ;
diagnostics.update.d  =d  ;
diagnostics.update.Lu =Lu ;
diagnostics.update.obstype=obstype;
diagnostics.update.bkgd=bkgd;
diagnostics.update.fcst=fcst;

% representer outputs: r.X_Y refers to sensitivity of model input variable Y
% to delta-perturbations of obs type X
diagnostics.r.H_h     =r_H(1:nx,:);
diagnostics.r.H_H0    =r_H(nx+1,:);
diagnostics.r.H_theta0=r_H(nx+2,:);
diagnostics.r.H_kadrag=r_H(nx+3,:);
diagnostics.r.v_h     =r_v(1:nx,:);
diagnostics.r.v_H0    =r_v(nx+1,:);
diagnostics.r.v_theta0=r_v(nx+2,:);
diagnostics.r.v_kadrag=r_v(nx+3,:);

% R-matrix outputs: R.X_Y refers to sensitivity of obs type Y to
% delta-perturbations of obs type X
diagnostics.R.H_H=R_HH;
diagnostics.R.H_v=R_Hv;
diagnostics.R.v_H=R_vH;
diagnostics.R.v_v=R_vv;

% additional diagnostic representers: rdiag.X_Y refers to sensitivity of
% model variable Y to delta-perturbations of obs type X
diagnostics.rdiag.H_H=r_HH;
diagnostics.rdiag.H_v=r_Hv;
diagnostics.rdiag.v_H=r_vH;
diagnostics.rdiag.v_v=r_vv;

% adjoint outputs: ad.X_Y refers to sensitivity of model input Y to
% delta-perturbations of obs type X
diagnostics.ad.H_h      =ad_H_h;
diagnostics.ad.H_H0     =ad_H_H0;
diagnostics.ad.H_theta0 =ad_H_theta0;
diagnostics.ad.H_ka_drag=ad_H_ka_drag;
diagnostics.ad.v_h      =ad_v_h;
diagnostics.ad.v_H0     =ad_v_H0;
diagnostics.ad.v_theta0 =ad_v_theta0;
diagnostics.ad.v_ka_drag=ad_v_ka_drag;