%   TV/TGV variational optical flow (stereo)
%
%   Author: Christoph Vogel
%
function [u,v,p] = solve_TV( I1, I2, A, u, v, p, par, constraints)

[M, N, C] = size(I1);
vv = cat( 1, u(:), v(:) );
y  = p(:);clear('p');

% add constraints at end
if exist('constraints', 'var' )
  numC = numel(constraints.uv(:))-1;
  y = cat(1, y, zeros(numel( constraints.uv), 1) );
end

normalizerR = 1./max(1, sum(abs(A),2));
normalizerC = 1./max(1, sum(abs(A),1)');

tt     = spdiags(normalizerR , 0, size(A,1), size(A,1));
tau    = reshape(normalizerC(1:N*M) , M,N);
ss     = spdiags(normalizerC, 0, size(A,2), size(A,2));
clear ('normalizerR');clear ('normalizerC');

ttA = tt*A;
ssA = ss*A';

if exist('constraints', 'var' )
  %note that after preconditioning new rhs is: ttrhs = tt*rhs(:);
  ttuv = tt* cat( 1, zeros(numel(y)-numC-1,1), constraints.uv(:));
  ttuv = ttuv(end-numC:end);% rest 0
end

clear ('tt');clear ('ss');clear ('A');

vv_=vv;
for w=1:par.warps
  
  u0 = u;
  v0 = v;
  
  % compute warping, etc.
  [I_x, I_y, I_t] = warping(I1, I2, u0, v0);
  I_grad_sqr = max(1e-09, I_x.^2 + I_y.^2);
  
  for i=1:par.maxits
    
    y = y + ttA*(vv_);

    if exist('constraints', 'var' )
      cons = y(end-numC:end) - ttuv;
      y1 = reshape(y(1:end-numC-1), [N*M, 4]);
      ys = sum(y1.^2, 2);
      consE = sqrt( ( cons(1:2:end).^2 + cons(2:2:end).^2));
      consE = repmat( max( 1, consE' ), [2,1]);
      cons = cons ./ consE(:);

      y = bsxfun(@times, y1, 1./max(1, sqrt(ys)));y=y(:);
      y = cat( 1, y, cons);
      clear('y1');
    else
      y1 = reshape(y, [N*M, 4]);
      ys = sum(y1.^2, 2);
      y = bsxfun(@times, y1, 1./max(1, sqrt(ys)));y=y(:);
      clear('y1');
    end
        
    vv_=vv;
    vv=vv-ssA*y;
    
    % extract u,v:
    u = reshape(vv(1:N*M), M,N);
    v = reshape(vv(1+N*M:2*N*M),M,N);
    
    % solve data term here == reprojection as well:    
    rho = I_t + (u-u0).*I_x + (v-v0).*I_y;
    idx = rho      < - par.lambda*tau.*I_grad_sqr;
    u(idx) = u(idx) + par.lambda*tau(idx).*I_x(idx);
    v(idx) = v(idx) + par.lambda*tau(idx).*I_y(idx);
    
    idx = rho      >   par.lambda*tau.*I_grad_sqr;
    u(idx) = u(idx) - par.lambda*tau(idx).*I_x(idx);
    v(idx) = v(idx) - par.lambda*tau(idx).*I_y(idx);
    
    idx = abs(rho) <=  par.lambda*tau.*I_grad_sqr;
    u(idx) = u(idx) - rho(idx).*I_x(idx)./I_grad_sqr(idx);
    v(idx) = v(idx) - rho(idx).*I_y(idx)./I_grad_sqr(idx);
    clear('idx');clear('rho');
    
    vv(1:N*M)       = u(:);
    vv(1+N*M:2*N*M) = v(:);
    
    vv_ = 2*vv-vv_;
    
  end
  % filter strong outliers
if ~exist('constraints', 'var' )  
  u = peakfilt(u);
  v = peakfilt(v);
end
  vv(1:N*M)       = u(:);
  vv(1+N*M:2*N*M) = v(:);
  vv_ = 2*vv-vv_;
end

if exist('constraints', 'var' )
  y(end-numC:end) = [];
end

p = reshape(y(1:N*M*4),M,N,4);