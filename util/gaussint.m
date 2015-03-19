% function [X,W]=gaussint(n,a,b,type,alfa,beta);
% The program returns the Gaussian integration sample points X and weights W. 
% INPUT  n         : number of sample points,  
%        type      : type of weight function and orthogonal polynomial and to be used, 
%        a,b       : integration interval limits, a < b,  
%        alfa,beta : powers for Jacobi weighting function (only for the type = 2) 
% OUTPUT X         : column vector of sample points (increasing order) 
%        W         : corresponding column vector of weights 
% 
% Type    Weight function     Orthogonal polynomial        Note 
%   1            1                 Legendre          
%   2   (x-a)^alfa*(b-x)^beta       Jacobi      alfa,beta > -1 & alfa^2+beta^2 > 0 
%   3          exp(-x)             Laguerre          [a,b] = [0,inf] 
%   4         exp(-x*x)             Hermite          [a,b] = [-inf,inf] 
% ---------------------------------------------------------------------------- 
% DEFAULT INPUT : gaussint(n)     = gaussint(n,-1,1,1) 
%                 gaussint(n,a,b) = gaussint(n,a,b,1) 
%                 gaussint(n,3)   = gaussint(n,0,inf,3) 
%                 gaussint(n,4)   = gaussint(n,-inf,inf,4) 
% ---------------------------------------------------------------------------- 
% 8 Feb 2003    : Jukka Sarvas, Electromagnetics Lab, Helsinki Univ. of Tech 
% ---------------------------------------------------------------------------- 
 
function [X,W]=gaussint(n,a,b,type,alfa,beta); 
if nargin <= 3 
   if nargin==2 
      type=a; 
   else 
      type=1; 
      if nargin==1 
         a=-1; 
         b=1; 
      end 
   end 
end 
m=1:n; 
if type==1 
   gamsq=(m-1).^2./((2*m-1).*(2*m-3)); 
   delta=zeros(1,n); 
   factor=2;   
   scalef=1; 
elseif type==2 
   if abs(1+alfa+beta)<1e-14 
      gamsq=[(m+beta-1).*(m+alfa-1)]./[(2*m+alfa+beta-2).^2]; 
      gamsq(1:2)=[0,2*gamsq(2)]; 
   else 
      gamsq=[4*(m-1).*(m+beta-1).*(m+alfa-1).*(m+alfa+beta-1)]./... 
         [(2*m+alfa+beta-1).*(2*m+alfa+beta-2).^2.*(2*m+alfa+beta-3)]; 
   end 
   delta=(alfa^2-beta^2)./((2*m+alfa+beta-2).*(2*m+alfa+beta)); 
   factor=2^(1+alfa+beta)*gamma(alfa+1)*gamma(beta+1)/gamma(alfa+beta+2); 
   scalef=((b-a)/2)^(alfa+beta); 
elseif type==3 
   gamsq=(m-1).^2;     
   delta=2*m-1; 
   factor=1; 
else 
   gamsq=(m-1)/2;   
   delta=zeros(1,n); 
   factor=sqrt(pi); 
end 
gam=sqrt(gamsq); 
J=diag(delta)+diag(gam(2:n),1)+diag(gam(2:n),-1); 
[V,X]=eig(J); 
X=diag(X);     
[X,I]=sort(X); 
W=zeros(n,1); 
for j=1:n 
   W(j)=factor*(V(1,j)/norm(V(:,j)))^2; 
end; 
W=W(I); 
if type < 3 
   W=.5*scalef*(b-a)*W; 
   X=.5*((b-a).*X+b+a); 
end 

