function [ad_h,ad_H0,ad_theta0,ad_ka_drag]=ad_waveModel(x,ad_H,ad_theta,ad_v,ad_k,bkgd)
%
% [ad_h,ad_H0,ad_theta0,ad_ka_drag]=ad_waveModel(x,ad_H,ad_theta,ad_v,ad_eps_r,ad_k,bkgd)
%
% AD-code for tl_waveModel.m.  Background state 'bkgd' can be a struct taken
% directly from output of waveModel.m, and 'tl_in' can be taken directly
% from output of tl_waveModel.m.
%
[g,alpha,beta,nu]=waveModelParams();

% grid
nx=length(x);
dx=diff(x(1:2));

% break out the bkgd vars
sigma=bkgd.sigma;
E    =bkgd.E    ;
Er   =bkgd.Er   ;
eps_b=bkgd.eps_b;
eps_r=bkgd.eps_r;
c    =bkgd.c    ;
cg   =bkgd.cg   ;
k    =bkgd.k    ;
h    =bkgd.h    ;
n    =bkgd.n    ;
theta=bkgd.theta;
sigma=bkgd.sigma;
x    =bkgd.x    ;
H    =bkgd.H    ;
dSxydx=bkgd.dSxydx;
Fy=bkgd.Fy;
v=bkgd.v;
H0=bkgd.H0;
theta0=bkgd.theta0;
gamma=bkgd.gamma;
Hm=bkgd.Hm;
Qb=bkgd.Qb;

ka_drag=bkgd.ka_drag;

refconst=sin(theta0)/c(nx);

%-----------------------------
% begin adjoint code
%-----------------------------

% initialize.  All these variables are assumed to never be directly
% observed, else they would be provided as inputs
zz=zeros(nx,1);
ad_c=zz;
ad_Er=zz;
ad_E=zz;
ad_eps_b=zz;
ad_eps_r=zz;
ad_cg=zz;
ad_h=zz;
ad_n=zz;
% ad_k=zz;
ad_Qb=zz;
ad_Hm=zz;
ad_refconst=0;
ad_gamma=0;

% bottom stress model following Ruessink et al. (2001), Feddersen et
% al. (2000)
a=1.16;  % empirical constant
Cd=0.015*(ka_drag./h).^(1/3);
urms=1.416*H*sigma./(4*sinh(k.*h));
% B=Cd./Fy.*v.^3./urms./sqrt(a^2+(v./urms).^2);  % old version of TL code, no mixing
B=a^2+(v./urms).^2;
dens = -urms.*Cd - v.^2.*Cd./urms./B;

% mixing operator
A=zeros(nx);
for i=2:nx-1
  A(i,i+[-1:1])=[1 -2 1]/dx^2*nu*h(i);
end
A(1,1:2)=[-2 1]/dx^2*nu*h(1);
A(nx,nx-1:nx)=[1 -2]/dx^2*nu*h(nx);

% % OLD: incorrect version for with-mixing case
% if(nu>0)
%   % v2: with mixing
%   %3b tl_v = inv(diag((1+B)./v)+A)*tl_N;
%   ad_N=inv(diag((1+B)./v)+A)'*ad_v;
%   ad_v=0;
%   %3a tl_N=tl_Fy./Fy - tl_Cd./Cd - tl_urms./urms.*(1-B);
%   ad_Fy=ad_N./Fy;
%   ad_Cd=-ad_N./Cd;
%   ad_urms=-ad_N./urms.*(1-B);
%   ad_N=0;
% NEW: corrected version for with-mixing case
if(nu>0)
  % v2: with mixing
  %3b tl_v = inv(diag(dens)+A)*tl_N;
  ad_N=inv(diag(dens)+A)'*ad_v;
else
  % v1: no mixing:
  %3b tl_v=tl_N./dens;
  ad_N=ad_v./dens;
end
ad_v=0;

%3a tl_N = tl_Fy./sqrt(B) + tl_urms.*(v.*Cd-v.^3.*Cd./(urms.^2.*B)) + tl_Cd.*v.*urms;
ad_Fy=ad_N./sqrt(B);
ad_urms=ad_N.*(v.*Cd-v.^3.*Cd./(urms.^2.*B));
ad_Cd=ad_N.*v.*urms;
ad_N=0;

%2 tl_urms=1.416*sigma*( tl_H./(4*sinh(k.*h)) ...
%                      -H./(4*sinh(k.*h).^2).*cosh(k.*h).*( tl_k.*h+k.*tl_h ) );
ad_H = ad_H + 1.416*sigma*( ad_urms./(4*sinh(k.*h)) );
ad_k = ad_k - 1.416*sigma*H./(4*sinh(k.*h).^2).*cosh(k.*h).*ad_urms.*h;
ad_h = ad_h - 1.416*sigma*H./(4*sinh(k.*h).^2).*cosh(k.*h).*ad_urms.*k;
ad_urms=0;

% note below: ka_drag is scalar, so sum over gridpoints.  Can see this by
% considering if the code was in a loop, then scalar ad_ka_drag would
% receive a contribution from each loop iteration
%1 tl_Cd=0.015*(1/3)*(ka_drag./h).^(-2/3).*(-ka_drag./h.^2.*tl_h+tl_ka_drag./h);
ad_h = ad_h + 0.015*(1/3)*(ka_drag./h).^(-2/3).*(-ka_drag./h.^2.*ad_Cd);
ad_ka_drag = sum(0.015*(1/3)*(ka_drag./h).^(-2/3).*(ad_Cd./h));  % if correcting ka_drag
ad_Cd=0;

% total force = radiation stress gradient + wind stress
% tl_Fy=tl_dSxydx;
ad_dSxydx=ad_Fy;
ad_Fy=0;

% radiation stress gradient
if(beta>0)
  % tl_dSxydx = -cos(theta)./c.*eps_r.*tl_theta ...
  %     +sin(theta)./c.^2.*eps_r.*tl_c ...
  %     -sin(theta)./c.*tl_eps_r;
  ad_theta=ad_theta-cos(theta)./c.*eps_r.*ad_dSxydx;
  ad_c=ad_c+sin(theta)./c.^2.*eps_r.*ad_dSxydx;
  ad_eps_r=ad_eps_r-sin(theta)./c.*ad_dSxydx;
  ad_dSxydx=0;
else
  % tl_dSxydx = -cos(theta)./c.*eps_b.*tl_theta ...
  %     +sin(theta)./c.^2.*eps_b.*tl_c ...
  %     -sin(theta)./c.*tl_eps_b;
  ad_theta=ad_theta-cos(theta)./c.*eps_b.*ad_dSxydx;
  ad_c=ad_c+sin(theta)./c.^2.*eps_b.*ad_dSxydx;
  ad_eps_b=ad_eps_b-sin(theta)./c.*ad_dSxydx;
  ad_dSxydx=0;
end

% stepping, explicit scheme
for i=1:(nx-1)

  % tl_H(i)=.5./sqrt(8/g*max(0,E(i)))*8/g.*tl_E(i);
  if(E(i)==0)
    ad_E(i)=0;
  else
    ad_E(i)=ad_E(i)+.5./sqrt(8/g*E(i))*8/g.*ad_H(i);
  end
  ad_H(i)=0;

  if(beta>0)

    % term 4
    nums1=2*Er(i+1)*c(i+1)*cos(theta(i+1));
    nums2=dx*(eps_b(i+1)-eps_r(i+1));
    denoms=2*c(i)*cos(theta(i));
    %4 tl_Er(i)=(tl_nums1+tl_nums2)/denoms ...
    %          - (nums1+nums2)/denoms^2*tl_denoms;
    ad_nums1=ad_Er(i)/denoms;
    ad_nums2=ad_Er(i)/denoms;
    ad_denoms=-(nums1+nums2)/denoms^2*ad_Er(i);
    ad_Er(i)=0;
    %3 tl_denoms=2*tl_c(i)*cos(theta(i)) ...
    %           - 2*c(i)*sin(theta(i))*tl_theta(i);
    ad_c(i)=ad_c(i)+2*cos(theta(i))*ad_denoms;
    ad_theta(i)=ad_theta(i)-2*c(i)*sin(theta(i))*ad_denoms;
    ad_denoms=0;
    %2 tl_nums2=dx*(tl_eps_b(i+1)-tl_eps_r(i+1));
    ad_eps_b(i+1)=ad_eps_b(i+1)+dx*ad_nums2;
    ad_eps_r(i+1)=ad_eps_r(i+1)-dx*ad_nums2;
    ad_nums2=0;
    %1 tl_nums1=2*tl_Er(i+1)*c(i+1)*cos(theta(i+1)) ...
    %          + 2*Er(i+1)*tl_c(i+1)*cos(theta(i+1)) ...
    %          - 2*Er(i+1)*c(i+1)*sin(theta(i+1))*tl_theta(i+1);
    ad_Er(i+1)=ad_Er(i+1)+2*c(i+1)*cos(theta(i+1))*ad_nums1;
    ad_c(i+1)=ad_c(i+1)+2*Er(i+1)*cos(theta(i+1))*ad_nums1;
    ad_theta(i+1)=ad_theta(i+1)-2*Er(i+1)*c(i+1)*sin(theta(i+1))*ad_nums1;
    ad_nums1=0;

    % term 3
    c1=2*g*sin(beta)/c(i+1);
    c2=-2*g*sin(beta)*Er(i+1)/c(i+1)^2;
    % tl_eps_r(i+1) = c1 * tl_Er(i+1) ...
    %     + c2 * tl_c(i+1);
    ad_Er(i+1)=ad_Er(i+1)+c1*ad_eps_r(i+1);
    ad_c(i+1) =ad_c(i+1) +c2*ad_eps_r(i+1);
    ad_eps_r(i+1)=0;

  end

  % term 2
  nums1=cg(i+1)*E(i+1)*cos(theta(i+1));
  nums2=eps_b(i+1)*dx;
  denoms=cg(i)*cos(theta(i));
  %4 tl_E(i) = tl_nums1/denoms ...
  %           - tl_nums2/denoms ...
  %           - (nums1-nums2)/denoms^2*tl_denoms;
  ad_nums1=ad_E(i)/denoms;  % note, consts initialized to zero
  ad_nums2=-ad_E(i)/denoms;
  ad_denoms=-(nums1-nums2)/denoms^2*ad_E(i);
  ad_E(i)=0;
  %3 tl_denoms=tl_cg(i)*cos(theta(i)) ...
  %           - cg(i)*sin(theta(i))*tl_theta(i);
  ad_cg(i)=ad_cg(i)+cos(theta(i))*ad_denoms;
  ad_theta(i)=ad_theta(i) - cg(i)*sin(theta(i))*ad_denoms;
  ad_denoms=0;
  %2 tl_nums2=tl_eps_b(i+1)*dx;
  ad_eps_b(i+1)=ad_eps_b(i+1)+ad_nums2*dx;
  ad_nums2=0;
  %1 tl_nums1=tl_cg(i+1)*E(i+1)*cos(theta(i+1)) ...
  %          + cg(i+1)*tl_E(i+1)*cos(theta(i+1)) ...
  %          - cg(i+1)*E(i+1)*sin(theta(i+1))*tl_theta(i+1);
  ad_cg(i+1)=ad_cg(i+1)+ad_nums1*E(i+1)*cos(theta(i+1));
  ad_E(i+1)=ad_E(i+1)+ad_nums1*cg(i+1)*cos(theta(i+1));
  ad_theta(i+1)=ad_theta(i+1)-ad_nums1*cg(i+1)*E(i+1)*sin(theta(i+1));
  ad_nums1=0;

  % term 1
  c1=alpha/4*g*(sigma/2/pi);
  % tl_eps_b(i+1)=c1*tl_Qb(i+1)*Hm(i+1)^2 ...
  %     + 2*c1*Qb(i+1)*Hm(i+1)*tl_Hm;
  ad_Qb(i+1)=ad_Qb(i+1)+c1*Hm(i+1)^2*ad_eps_b(i+1);
  ad_Hm(i+1)=ad_Hm(i+1)+2*c1*Qb(i+1)*Hm(i+1)*ad_eps_b(i+1);
  ad_eps_b(i+1)=0;

  % fraction of breaking waves, non-implicit approximation from SWAN code
  B=H(i+1)/Hm(i+1);
  if(B<=.5)
    Qo=0;
  else
    Qo=(2*B-1)^2;
  end
  if(.2<B<=1)
    args=(Qo-1)/B^2;
    nums=Qo-exp(args);
    dens=B^2-exp(args);
    % Qb(i+1)=Qo-B^2*nums/dens;
    % tl_Qb(i+1)=tl_Qo ...
    %     - 2*B*nums/dens*tl_B ...
    %     - B^2*tl_nums/dens ...
    %     + B^2*nums/dens^2*tl_dens;
    % note: adjoint temporary scalar vars not initialized, assume zero to start
    ad_Qo=ad_Qb(i+1);
    ad_B=-ad_Qb(i+1)*2*B*nums/dens;
    ad_nums=-B^2*ad_Qb(i+1)/dens;
    ad_dens=B^2*nums/dens^2*ad_Qb(i+1);
    ad_Qb(i+1)=0;
    % tl_dens=2*B*tl_B - exp(args)*tl_args;
    ad_B=ad_B+2*B*ad_dens;
    ad_args=-exp(args)*ad_dens;
    ad_dens=0;
    % tl_nums=tl_Qo - exp(args)*tl_args;
    ad_Qo=ad_Qo+ad_nums;
    ad_args=ad_args-exp(args)*ad_nums;
    ad_nums=0;
    % tl_args=tl_Qo/B^2-2*(Qo-1)/B^3*tl_B;
    ad_Qo=ad_Qo+ad_args/B^2;
    ad_B=ad_B-2*(Qo-1)/B^3*ad_args;
    ad_args=0;
  else
    % tl_Qb(i+1)=0;
    ad_Qb(i+1)=0;
  end

  if(B<=.5)
    % tl_Qo=0;
    ad_Qo=0;
  else
    % tl_Qo=2*(2*B-1).*2*tl_B;
    ad_B=ad_B+2*(2*B-1).*2*ad_Qo;
    ad_Qo=0;
  end

  % tl_B=tl_H(i+1)/Hm(i+1)-H(i+1)/Hm(i+1)^2*tl_Hm(i+1);
  ad_H(i+1)=ad_H(i+1)+ad_B/Hm(i+1);
  ad_Hm(i+1)=ad_Hm(i+1)-H(i+1)/Hm(i+1)^2*ad_B;
  ad_B=0;

  % max wave height
  % Hm=0.88./k.*tanh(gamma/0.88.*k.*h);
  tharg=gamma/0.88.*k(i+1).*h(i+1);
  % tl_Hm(i+1)=0.88*( -1./k(i+1).^2.*tanh(tharg).*tl_k(i+1) ...
  %              + 1./k(i+1).*sech(tharg).^2.*tl_tharg );
  ad_k(i+1)=ad_k(i+1)-0.88*ad_Hm(i+1)./k(i+1).^2.*tanh(tharg);
  ad_tharg=ad_Hm(i+1)*0.88./k(i+1).*sech(tharg).^2;  % assume init = 0
  ad_Hm(i+1)=0;
  % tl_tharg=gamma/0.88*( tl_k(i+1).*h(i+1) + k(i+1)*tl_h(i+1) ) ...
  %          + tl_gamma/0.88.*k(i+1).*h(i+1);
  ad_k(i+1)=ad_k(i+1)+tharg./k(i+1)*ad_tharg;
  ad_h(i+1)=ad_h(i+1)+tharg./h(i+1)*ad_tharg;
  ad_gamma=ad_gamma+1/0.88.*k(i+1).*h(i+1)*ad_tharg;
  ad_tharg=0;

  % refraction
  const=1./sqrt(1-(c(i).*refconst).^2);
  % tl_theta(i)=const.*( refconst.*tl_c(i) + tl_refconst.*c(i) );
  ad_c(i)=ad_c(i)+const.*refconst.*ad_theta(i);
  ad_refconst=ad_refconst+const.*c(i)*ad_theta(i);
  ad_theta(i)=0;

end
% tl_theta(nx)=1./sqrt(1-(c(nx).*refconst).^2).*(refconst.*tl_c(nx) + tl_refconst.*c(nx) );
ad_c(nx)=ad_c(nx)+1./sqrt(1-(c(nx).*refconst).^2).*refconst.*ad_theta(nx);
ad_refconst=ad_refconst+1./sqrt(1-(c(nx).*refconst).^2).*ad_theta(nx).*c(nx);
ad_theta(nx)=0;
% tl_H(nx)=tl_H0;
ad_H0=ad_H(nx);
% tl_E(nx)=g/8*2*H0*tl_H0;
ad_H0=ad_H0+g/8*2*H0*ad_E(nx);
% tl_Er(nx)=0;

% gamma calculated based on deep water wave steepness (s0) following Battjes
% and Stive (1985), and also used at FRF by Ruessink et al. (2001)
L0=g/(2*pi*(sigma/2/pi)^2);
s0=H0/L0;
gamma=0.5+0.4*tanh(33*s0);
%2 tl_gamma=0.4*sech(33*s0).^2.*33*tl_s0;
ad_s0=0.4*sech(33*s0).^2.*33*ad_gamma;
ad_gamma=0;
%1 tl_s0=tl_H0/L0;
ad_H0=ad_H0+ad_s0/L0;
ad_s0=0;

% dispersion

%5 tl_refconst=cos(theta0)/c(nx)*tl_theta0-sin(theta0)/c(nx)^2*tl_c(nx);
ad_theta0=cos(theta0)/c(nx)*ad_refconst;
ad_c(nx)=ad_c(nx)-sin(theta0)/c(nx)^2*ad_refconst;
ad_refconst=0;

%4 tl_cg=tl_n.*c+n.*tl_c;
ad_n=ad_n+c.*ad_cg;
ad_c=ad_c+n.*ad_cg;
ad_cg=0;

%3 tl_n = tl_k.*h./sinh(2*k.*h) + k.*tl_h./sinh(2*k.*h) ...
%        - k.*h./sinh(2*k.*h).^2.*cosh(2*k.*h)*2.*(tl_k.*h+k.*tl_h);
ad_k=ad_k+h./sinh(2*k.*h).*ad_n;
ad_h=ad_h+k./sinh(2*k.*h).*ad_n;
ad_k=ad_k-k.*h./sinh(2*k.*h).^2.*cosh(2*k.*h)*2.*h.*ad_n;
ad_h=ad_h-k.*h./sinh(2*k.*h).^2.*cosh(2*k.*h)*2.*k.*ad_n;
ad_n=0;

%2 tl_c=-sigma./k.^2.*tl_k;
ad_k=ad_k-sigma./k.^2.*ad_c;
ad_c=0;

%1 tl_k=-tl_h.*k.^2.*sech(k.*h).^2./(tanh(k.*h)+k.*h.*sech(k.*h).^2);
ad_h=ad_h-k.^2.*sech(k.*h).^2./(tanh(k.*h)+k.*h.*sech(k.*h).^2).*ad_k;
ad_k=0;

ad_h=ad_h(:);