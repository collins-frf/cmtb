function out=waveModel(x,H0,theta0,out)
%
% out=waveModel(x,H0,theta0,in)
%
% Wave energy balance equation solver, explicit spatial stepping scheme.
%
% Breaking dissipation (eps_b) using TG1983
% Roller energy (Er) and roller dissipation (eps_r) following Reniers & Battjes (1996)
%
% NOTE: input theta0 in radians
%
% 'in' will be appended/overwritten with new variables to create 'out', and
% should contain the following minimal inputs:
%
% in.{h,H0,theta0,sigma,tauw}
%
% note, tauw is optional, alongshore component of wind stress in m2/s2 units
%

h     =out.h     ;
sigma =out.sigma ;
ka_drag=out.ka_drag;

% wind stress implemented later, so use as optional argument to ensure
% backwards compatibility
if(isfield(out,'tauw'))
  tauw=out.tauw;
else
  tauw=0;
end

[g,alpha,beta,nu,gammaType]=waveModelParams();

% grid
nx=length(x);
dx=diff(x(1:2));

% dispersion
% k=fsolve(@(k)sigma^2-g*k.*tanh(k.*h),sigma./sqrt(g*h),optimset('Display','off'));
k=nan*h;
for i=1:nx
  k(i)=fzero(@(k)sigma^2-g*k.*tanh(k.*h(i)),sigma./sqrt(g*h(i)),optimset('Display','off'));
end
c=max(0,real(sigma./k));
n=.5*(1+2*k.*h./sinh(2*k.*h));
cg=n.*c;
refconst=sin(theta0)/c(nx);

% gamma can be either calculated based on deep water wave steepness (s0)
% following Battjes and Stive (1985) (also used by Ruessink et al., 2001),
% or based on the empirical fit obtained for duck94 by Ruessink et
% al. (2003).
if(gammaType==2001)
  L0=g/(2*pi*(omega/2/pi)^2);
  s0=H0/L0;
  gamma=0.5+0.4*tanh(33*s0);
  gamma=ones(nx,1)*gamma;
elseif(gammaType==2003)
  gamma=0.76*k.*h+0.29;
end

% refraction
theta=asin(c.*refconst);

% stepping, explicit scheme
E=zeros(nx,1);
Er=zeros(nx,1);
eps_b=zeros(nx,1);
eps_r=zeros(nx,1);
E(nx)=g/8*H0^2;
Er(nx)=0;
H(nx)=H0;
theta(nx)=theta0; %asin(c(nx).*refconst);
for i=(nx-1):-1:1

  % max wave height
  tharg=gamma(i+1)/0.88.*k(i+1).*h(i+1);
  Hm(i+1)=0.88./k(i+1).*tanh(tharg);

  % fraction of breaking waves, non-implicit approximation from SWAN code
  B=H(i+1)/Hm(i+1);
  if(B<=.5)
    Qo=0;
  else
    Qo=(2*B-1)^2;
  end
  if(B<=.2)
    Qb(i+1)=0;
  elseif(.2<B<=1)
    args=(Qo-1)/B^2;
    nums=Qo-exp(args);
    dens=B^2-exp(args);
    Qb(i+1)=Qo-B^2*nums/dens;
  else
    Qb(i+1)=1;
  end
  c1=alpha/4*g*(sigma/2/pi);
  eps_b(i+1)=c1*Qb(i+1)*Hm(i+1)^2;

  nums1=cg(i+1)*E(i+1)*cos(theta(i+1));
  nums2=eps_b(i+1)*dx;
  denoms=cg(i)*cos(theta(i));
  E(i)=(nums1-nums2)/denoms;

  if(beta>0)
    eps_r(i+1)=2*g*Er(i+1)*sin(beta)/c(i+1);
    Er(i)=(2*Er(i+1)*c(i+1)*cos(theta(i+1))+dx*(eps_b(i+1)-eps_r(i+1)))/(2*c(i)*cos(theta(i)));
    if(Er(i)<0)
      Er(i)=0;
    end
  end

  if(E(i)<.001)
    E(i)=.001;
  end
  H(i)=sqrt(8/g*E(i));

end
c=c(:);
cg=cg(:);
k=k(:);
n=n(:);
theta=theta(:);
H=H(:);

% radiation stress gradient
if(beta>0)  % roller
  dSxydx = -sin(theta)./c.*eps_r;
else
  dSxydx=-sin(theta)./c.*eps_b;
end
dSxydx(dSxydx==0)=1e-6;  % avoid singularity in TL model

% bottom stress model following Ruessink et al. (2001), Feddersen et al. (2000)

% total force = radiation stress gradient + wind stress
Fy=dSxydx+tauw;

% v1: analytical solution, no mixing
a=1.16;  % empirical constant
Cd=0.015*(ka_drag./h).^(1/3);
urms=1.416*H.*sigma./(4*sinh(k.*h));
v2 = sqrt( (a*Cd.*urms).^4 + 4*(Cd.*Fy).^2 )./(2*Cd.^2) - (a*urms).^2/2;
v=sqrt(v2).*sign(-Fy);

% mixing operator
A=zeros(nx);
for i=2:nx-1
  A(i,i+[-1:1])=[1 -2 1]/dx^2*nu*h(i);
end
A(1,1:2)=[-2 1]/dx^2*nu*h(1);
A(nx,nx-1:nx)=0;  % [1 -2]/dx^2*nu*h(nx);

% v2: nonlinear solution with mixing.
v0=v;
v = fsolve(@(v)Fy + Cd.*urms.*v.*sqrt(a^2+(v./urms).^2) - A*v,v0,optimset('Display','off'));

% outputs struct
out.E    =E    ;
out.Er   =Er   ;
out.eps_b=eps_b;
out.eps_r=eps_r;
out.c=c;
out.cg=cg;
out.k=k;
out.h=h;
out.n=n;
out.theta=theta;
out.sigma=sigma;
out.x=x;
out.H=H;
out.gamma=gamma;
out.Hm=Hm;
out.Qb=Qb;
out.dSxydx=dSxydx;
out.Fy=Fy;
out.v=real(v);
out.ka_drag=ka_drag;
