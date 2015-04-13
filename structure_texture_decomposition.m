%
%   ROF/structure texture decomposition
%
function outI = structure_texture_decomposition(I, par, maxIterations)

factor  = par.structure_texture_factor;

if ~exist('maxIterations','var')
    maxIterations = 20;
end

backupI = I;
outI    = I;

for iIm = 1:size(I,3)
  
    I = squeeze(backupI(:,:,iIm));
  
    [M,N] = size(I);
    pu1 = zeros ( M, N );
    pu2 = zeros ( M, N );
    sigma = 0.5;
    tau   = 0.25;
    lambda = 4.0;%2.0; % lambda = 'strength' of data term
    
    u_ = I;
    u  = u_;  
    % structure to texture decomposition by ROF model
    
    for i=1:maxIterations

      % compute derivatives
      u_x = dxp(u_);
      u_y = dyp(u_);

      u_ = u;% u_ = u_n
      
      % update dual variable
      pu1 = pu1 + sigma*u_x;
      pu2 = pu2 + sigma*u_y;

      % reprojection to |pu| <= 1
      reprojection = max(1.0, sqrt(pu1.^2 + pu2.^2));
      pu1 = pu1./reprojection;
      pu2 = pu2./reprojection;
      
      % compute divergence
      div_p = dxm(pu1) + dym(pu2);

      % min_v prox_tau,D(v) = 1/tau (v-u)^2 + 1/lambda*(v-f)^2 
      % 2(v-f) + 2/tau(v-u) = 0 <-> 2*(1+1/tau) v = 2*(f+1/tau u)
      % v = (f+1/tau*u)  / (1+1/tau)
      % thus 
      % 1/(tau)( (2f+2u/tau)/(2+2/tau) -u)^2 + ((2f+2u/tau)  / (2+2/tau)-f)^2 
      % which is
      % u = (u+2tau f)/(1+2 tau)
      
      % compute u_n+1
      u = (u + tau*div_p + (lambda*2*tau) * I)./(1+lambda*2*tau);
      
      % extapolation:
      theta = 1/sqrt(1+4*tau); tau = theta*tau; sigma = sigma/theta;
      u_ = u+theta*(u-u_);
      
    end
    
   	s1 = mean(mean(I));
    I = (1-factor) * u + (I - u) * factor;
  
    s2 = mean(mean(I));
    I = I * s1/s2;

    I = scale_image(I, 0.0, 1.0);
  
  outI(:,:,iIm) = I;
  
%  outI(:,:,iIm) = u;
  
end