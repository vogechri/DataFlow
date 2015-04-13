% structure texture decomposition according to 
% Structure-Texture Image Decomposition--Modeling, Algorithms, and Parameter Selection IJCV 2006
% Author: Christoph Vogel
function I = structure_texture_decompositionMask(I, factor, mask)

if factor == 0.5
  return;
end

  [M,N] = size(I);
  maxIterations = 50;
  % the higher theta the more it is a pure approximation to I
  theta = 0.05;

  u  = zeros(size(I(:)));

if exist('mask','var')
  A = System_Matrix_TV( 1, M,N, 1, mask );
else
  A = System_Matrix_TV( 1, M,N, 1 );
end
  normalizerR = 1./max(1, sum(abs(A),2));
  normalizerC = 1./max(1, sum(abs(A),1)');

  tt     = spdiags(normalizerR , 0, size(A,1), size(A,1));
  tau    = reshape(normalizerC(1:N*M) , M,N);
  ss     = spdiags(normalizerC, 0, size(A,2), size(A,2));
  clear ('normalizerR');clear ('normalizerC');

  ttA = tt*A;
  ssA = ss*A';
  clear ('tt');clear ('ss');clear ('A');
  
  y  = zeros(2*numel(u), 1);
  u_ = u;
  for i=1:maxIterations
  
    y = y + ttA*(u_(:));
    y1 = reshape(y, [N*M, 2]);
    ys = sum(y1.^2, 2);
    y = bsxfun(@times, y1, 1./max(1, sqrt(ys)));y=y(:);
    
    u_= u;
    u = u-ssA*y;
 
    % compute u
    u = (theta.*I(:) + tau(:).*u)./(theta+tau(:));
    u_ = 2*u-u_;
  end

  u=reshape(u,size(I));
  I = (1-factor) * I + (I - u) * factor;
  %  I = scale_image(I, 0.0, 1.0);  
end