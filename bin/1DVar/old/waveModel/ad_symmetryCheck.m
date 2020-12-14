%
% tests adjoint code: F'*AD*TL*F should be symmetric, +'ve definite
%
clear

% get NL solution (background)
x=[100:10:1000]';
bkgd.h=.01*x;
bkgd.H0=1;
bkgd.theta0=deg2rad(10);
bkgd.sigma=2*pi/10;
bkgd.ka_drag=0.015;
nx=length(x);
bkgd=waveModel(x,bkgd.H0,bkgd.theta0,bkgd);

% apply TL and ADJ models for n instances of random forcing F (actually, IC
% perturbations)
eps = 0.05;
n=10;
F = eps*rand(nx,n);
F(end,:)=F(1,:);
for i=1:n
  % TL model: u=TL*F
  [tl_H,tl_theta,tl_v,tl_k]=tl_waveModel(x,F(:,i),0,0,0,bkgd);
  % ADJ model: g=ADJ*(TL*F)
  g(:,i)=ad_waveModel(x,tl_H,tl_theta,tl_v,tl_k,bkgd);
end

% test whether result makes sense
A = F'*g;

(A-A')/sum(diag(A))  % should be zeros to within roundoff

%pcolor(A)
format short g
eig(A)  % should be +'ve
