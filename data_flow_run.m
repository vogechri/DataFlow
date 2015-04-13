%%%   TV/TGV variational optical flow
%%%
%%% An Evaluation of Data Costs for Optical Flow, 
%%% C. Vogel, S. Roth, K. Schindler, 
%%% In: German Conference on Pattern Recognition (GCPR), Saarbrücken, Germany, 2013

function [u, v, w, p, q] = data_flow_run(I1, I2, u, v, w, p, q, par, mask, level_constraints )

% TGV only parameters, relative weight is important here, so fix alpha1 =1
alpha0 = 2; % 1.5 parameter defining the influence of second integral for TGV
alpha1 = 1; % 1   parameter defining the influence of first integral for TGV

[M, N, C] = size(I1);

if ~exist('level_constraints', 'var')
  level_constraints.p = [];
  level_constraints.uv = [];
else
  level_constraints.uv = bsxfun(@times, level_constraints.uv', level_constraints.weights);
  % remove points too close at boundary:
  ids = (level_constraints.p(:,1)<2 | level_constraints.p(:,2)<2 | level_constraints.p(:,1)>N-1 | level_constraints.p(:,2)>M-1);
  
  level_constraints.uv = level_constraints.uv(:,~ids);
  level_constraints.p = level_constraints.p(~ids,:);
  level_constraints.weights = level_constraints.weights(~ids);
  
  fu=u( sub2ind( [M,N], round( level_constraints.p(:,2)), round( level_constraints.p(:,1)) )) + level_constraints.p(:,1);
  fv=v( sub2ind( [M,N], round( level_constraints.p(:,2)), round( level_constraints.p(:,1)) )) + level_constraints.p(:,2);
  
  ids = (fu<2 | fv<2 | fu>N-1 | fv>M-1);
  level_constraints.uv = level_constraints.uv(:,~ids);
  level_constraints.p = level_constraints.p(~ids,:);
  level_constraints.weights = level_constraints.weights(~ids);  
end

if ~exist('mask', 'var')
  mask = false(size(I1));
else
  u(mask) = 0;
  v(mask) = 0;
  w(repmat(mask, [1,1,size(w,3)])) = 0;
  p(repmat(mask, [1,1,size(p,3)])) = 0;
  q(repmat(mask, [1,1,size(q,3)])) = 0;
end

  if par.doTV
    if par.edgeFilter
      g = edge_Tensor(I1);
      [yids, xids, vids] = Matrix_data_TV( alpha1, M,N, par, mask, level_constraints, g );
    else
      [yids, xids, vids] = Matrix_data_TV( alpha1, M,N, par, mask, level_constraints );
    end
  else % TGV:
    if par.edgeFilter
      g = edge_Tensor(I1); 
      [yids, xids, vids] = Matrix_data_TGV( alpha1, alpha0, M,N, par, mask, level_constraints, g );      
    else
      [yids, xids, vids] = Matrix_data_TGV( alpha1, alpha0, M,N, par, mask, level_constraints );
    end
  end
  
  if ~isempty(xids)
%    G = sparse(yids, xids, vids);
    
    [u,v,w,p,q] = ofmex( I1,I2, par.ring, u,v, par.lambda, yids, xids, vids, ...
                         max(yids), max(xids), cat(1, u(:), v(:), w(:) ), ...
                         cat(1, p(:), q(:)), par.warps, par.maxits, par.cEps, ...
                         par.doStereo, par.dataTerm, par.doTV, ...
                         level_constraints.uv(:), par.medFilt );
  end

if exist('mask', 'var')
  u(mask) = 0;
  v(mask) = 0;
  w(repmat(mask, [1,1,size(w,3)])) = 0;
  p(repmat(mask, [1,1,size(p,3)])) = 0;
  q(repmat(mask, [1,1,size(q,3)])) = 0;
end
