function [tl_H,tl_theta,tl_v,tl_k]=tl_waveModel(x,tl_h,tl_H0,tl_theta0,tl_ka_drag,bkgd)
%
% [tl_H,tl_theta,tl_v,tl_k]=tl_waveModel(x,tl_h,tl_H0,tl_theta0,tl_ka_drag,bkgd)
%
% TL-code for waveModel.m.  Background state 'bkgd' can be a struct taken
% directly from output of waveModel.m
%

[g,alpha,beta,nu]=waveModelParams();

% grid
nx=length(x);
dx=diff(x(1:2));

% break out the bkgd vars
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

%-----------------------------
% begin tangent-linear code
%-----------------------------

% dispersion
tl_k=-tl_h.*k.^2.*sech(k.*h).^2./(tanh(k.*h)+k.*h.*sech(k.*h).^2);
tl_c=-sigma./k.^2.*tl_k;
% n= .5* + k.*h./sinh(2*k.*h);
tl_n = tl_k.*h./sinh(2*k.*h) + k.*tl_h./sinh(2*k.*h) ...
       - k.*h./sinh(2*k.*h).^2.*cosh(2*k.*h)*2.*(tl_k.*h+k.*tl_h);
tl_cg=tl_n.*c+n.*tl_c;
refconst=sin(theta0)/c(nx);
tl_refconst=cos(theta0)/c(nx)*tl_theta0-sin(theta0)/c(nx)^2*tl_c(nx);

% gamma calculated based on deep water wave steepness (s0) following Battjes
% and Stive (1985), and also used at FRF by Ruessink et al. (2001)
L0=g/(2*pi*(sigma/2/pi)^2);
s0=H0/L0;
tl_s0=tl_H0/L0;
gamma=0.5+0.4*tanh(33*s0);
tl_gamma=0.4*sech(33*s0).^2.*33*tl_s0;

% stepping, explicit scheme
tl_E=zeros(nx,1);
tl_Er=zeros(nx,1);
tl_eps_b=zeros(nx,1);
tl_eps_r=zeros(nx,1);
tl_E(nx)=g/8*2*H0*tl_H0;
tl_Er(nx)=0;
tl_H(nx)=tl_H0;
tl_theta(nx)=tl_theta0;
for i=(nx-1):-1:1

  % refraction
  % theta(i)=asin(c(i).*refconst);
  tl_theta(i)=1./sqrt(1-(c(i).*refconst).^2).*( refconst.*tl_c(i) + tl_refconst.*c(i) );

  % max wave height
  tharg=gamma/0.88.*k(i+1).*h(i+1);
  tl_tharg=gamma/0.88*( tl_k(i+1).*h(i+1) + k(i+1)*tl_h(i+1) ) ...
           + tl_gamma/0.88.*k(i+1).*h(i+1);
  % Hm(i+1)=0.88./k(i+1).*tanh(tharg);
  tl_Hm(i+1)=0.88*( -1./k(i+1).^2.*tanh(tharg).*tl_k(i+1) ...
               + 1./k(i+1).*sech(tharg).^2.*tl_tharg );

  % fraction of breaking waves, non-implicit approximation from SWAN code
  B=H(i+1)/Hm(i+1);
  tl_B=tl_H(i+1)/Hm(i+1)-H(i+1)/Hm(i+1)^2*tl_Hm(i+1);
  if(B<=.5)
    Qo=0;
    tl_Qo=0;
  else
    Qo=(2*B-1)^2;
    tl_Qo=2*(2*B-1).*2*tl_B;
  end
  if(B<=.2)
    tl_Qb(i+1)=0;
  elseif(.2<B<=1)
    args=(Qo-1)/B^2;
    tl_args=tl_Qo/B^2-2*(Qo-1)/B^3*tl_B;
    nums=Qo-exp(args);
    tl_nums=tl_Qo-exp(args)*tl_args;
    dens=B^2-exp(args);
    tl_dens=2*B*tl_B-exp(args)*tl_args;
    Qb(i+1)=Qo-B^2*nums/dens;
    tl_Qb(i+1)=tl_Qo-2*B*tl_B*nums/dens ...
        -B^2*tl_nums/dens + B^2*nums/dens^2*tl_dens;
  else
    tl_Qb(i+1)=0;
  end

  % term 1
  c1=alpha/4*g*(sigma/2/pi);
  % eps_b(i+1)=c1*Qb(i+1)*Hm(i+1)^2;
  tl_eps_b(i+1)=c1*tl_Qb(i+1)*Hm(i+1)^2 ...
      + 2*c1*Qb(i+1)*Hm(i+1)*tl_Hm(i+1);

  % term 2
  nums1=cg(i+1)*E(i+1)*cos(theta(i+1));
  tl_nums1=tl_cg(i+1)*E(i+1)*cos(theta(i+1)) ...
           + cg(i+1)*tl_E(i+1)*cos(theta(i+1)) ...
           - cg(i+1)*E(i+1)*sin(theta(i+1))*tl_theta(i+1);
  nums2=eps_b(i+1)*dx;
  tl_nums2=tl_eps_b(i+1)*dx;
  denoms=cg(i)*cos(theta(i));
  tl_denoms=tl_cg(i)*cos(theta(i)) ...
            - cg(i)*sin(theta(i))*tl_theta(i);
  % E(i)=(nums1-nums2)/denoms;
  tl_E(i) = (tl_nums1-tl_nums2)/denoms ...
            - (nums1-nums2)/denoms^2*tl_denoms;

  if(beta>0)

    % term 3
    % eps_r(i+1)=2*g*Er(i+1)*sin(beta)/c(i+1);
    tl_eps_r(i+1)=2*g*tl_Er(i+1)*sin(beta)/c(i+1) ...
        - 2*g*Er(i+1)*sin(beta)/c(i+1)^2*tl_c(i+1);

    % term 4
    nums1=2*Er(i+1)*c(i+1)*cos(theta(i+1));
    nums2=dx*(eps_b(i+1)-eps_r(i+1));
    denoms=2*c(i)*cos(theta(i));
    tl_nums1=2*tl_Er(i+1)*c(i+1)*cos(theta(i+1)) ...
             + 2*Er(i+1)*tl_c(i+1)*cos(theta(i+1)) ...
             - 2*Er(i+1)*c(i+1)*sin(theta(i+1))*tl_theta(i+1);
    tl_nums2=dx*(tl_eps_b(i+1)-tl_eps_r(i+1));
    tl_denoms=2*tl_c(i)*cos(theta(i)) ...
              - 2*c(i)*sin(theta(i))*tl_theta(i);
    % Er(i)=(nums1+nums2)/denoms;
    tl_Er(i)=(tl_nums1+tl_nums2)/denoms ...
             - (nums1+nums2)/denoms^2*tl_denoms;

  end

  % H(i)=sqrt(8/g*E(i));
  if(E(i)==0)
    tl_H(i)=0;
  else
    tl_H(i)=.5./sqrt(8/g*E(i))*8/g.*tl_E(i);
  end

end
tl_c=tl_c(:);
tl_H=tl_H(:);
tl_theta=tl_theta(:);

% radiation stress gradient
if(beta>0)
  % dSxydx = -sin(theta)./c.*eps_r;
  tl_dSxydx = -cos(theta)./c.*eps_r.*tl_theta ...
      +sin(theta)./c.^2.*eps_r.*tl_c ...
      -sin(theta)./c.*tl_eps_r;
else
  % dSxydx = -sin(theta)./c.*eps_b;
  tl_dSxydx = -cos(theta)./c.*eps_b.*tl_theta ...
      +sin(theta)./c.^2.*eps_b.*tl_c ...
      -sin(theta)./c.*tl_eps_b;
end

% total force = radiation stress gradient + wind stress
% Fy=dSxydx+tauw;
tl_Fy=tl_dSxydx;

% define mixing operator
A=zeros(nx);
for i=2:nx-1
  A(i,i+[-1:1])=[1 -2 1]/dx^2*nu*h(i);
end
A(1,1:2)=[-2 1]/dx^2*nu*h(1);
A(nx,nx-1:nx)=[1 -2]/dx^2*nu*h(nx);

% bottom stress model following Ruessink et al. (2001), Feddersen et
% al. (2000).  To get TL model, differentiate the eqn for v (i.e., the
% fsolve() line in waveModel.m) on both sides, then solve for
% tl_v
a=1.16;  % empirical constant
Cd=0.015*(ka_drag./h).^(1/3);
tl_Cd=0.015*(1/3)*(ka_drag./h).^(-2/3).*(-ka_drag./h.^2.*tl_h+tl_ka_drag./h);
urms=1.416*H*sigma./(4*sinh(k.*h));
tl_urms=1.416*sigma*( tl_H./(4*sinh(k.*h)) ...
                      -H./(4*sinh(k.*h).^2).*cosh(k.*h).*( tl_k.*h+k.*tl_h ) );
B=a^2+(v./urms).^2;
tl_N = tl_Fy./sqrt(B) + tl_urms.*(v.*Cd-v.^3.*Cd./(urms.^2.*B)) + tl_Cd.*v.*urms;
dens = -urms.*Cd - v.^2.*Cd./urms./B;
if(nu==0)
  tl_v=tl_N./dens;
else
  tl_v = inv(diag(dens)+A)*tl_N;  % tl_N = (dens + A) * tl_v
end





return;
%-----------------------------------
% OLD: incorrect derivation of tl_v, did not correctly incorporate mixing
%-----------------------------------

% bottom stress model following Ruessink et al. (2001), Feddersen et
% al. (2000).  To get TL model, differentiate the eqn for tau_b=Fy on both
% sides, then solve for tl_v... double-checked algebra on this, and verified
% that tl_v is consistent with a perturbed nonlinear model (but this check
% only applied to the version without mixing)
a=1.16;  % empirical constant
Cd=0.015*(ka_drag./h).^(1/3);
tl_Cd=0.015*(1/3)*(ka_drag./h).^(-2/3).*(-ka_drag./h.^2.*tl_h+tl_ka_drag./h);
urms=1.416*H*sigma./(4*sinh(k.*h));
tl_urms=1.416*sigma*( tl_H./(4*sinh(k.*h)) ...
                      -H./(4*sinh(k.*h).^2).*cosh(k.*h).*( tl_k.*h+k.*tl_h ) );
B=Cd./Fy.*v.^3./urms./sqrt(a^2+(v./urms).^2);
nums = tl_Fy./Fy ...
       - tl_Cd./Cd ...
       - tl_urms./urms.*(1-B);
dens = (1+B)./v;
tl_v = nums./dens;

% v2: with mixing operator.  Note previously without mixing I had worked out
% the TL form of tau_b(v)=Fy, to get
%
% nums = dens.*tl_v
%      = diag(dens)*tl_v.
%
% Now just add in the mixing operator A to the same derivation, to get
%
% nums = ( diag(dens) + A )*tl_v.
% ---> tl_v = inv(...)*nums
%
A=zeros(nx);
for i=2:nx-1
  A(i,i+[-1:1])=[1 -2 1]/dx^2*nu*h(i);
end
A(1,1:2)=[-2 1]/dx^2*nu*h(1);
A(nx,nx-1:nx)=[1 -2]/dx^2*nu*h(nx);
if(nu==0)
  tl_v=nums./dens;
else
  keyboard;
  tl_v = pinv(diag(dens)+A)*nums;
end
