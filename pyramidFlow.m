%%%   TV/TGV variational optical flow
%%%
%%% An Evaluation of Data Costs for Optical Flow, 
%%% C. Vogel, S. Roth, K. Schindler, 
%%% In: German Conference on Pattern Recognition (GCPR), Saarbrücken, Germany, 2013
%

function flow_out = pyramidFlow ( I1, I2, p, startLev, constraints )

[M, N, C] = size(I1);
I1 = scale_image(double(I1), 0, 1);
I2 = scale_image(double(I2), 0, 1);

if ~exist('startLev', 'var')
  startLev=50;% the minimal size for an image in the pyramid
end

width_Pyramid  = cell(1,1);
height_Pyramid = cell(1,1);

%N*p.pyramid_factor^levelX = 28;
levelX = log(startLev/N)/log(p.pyramid_factor);
levelY = log(startLev/M)/log(p.pyramid_factor);
levels = round(min(levelY, levelX));
pyramid_factor_x = exp(log(startLev/N)/levels);
pyramid_factor_y = exp(log(startLev/M)/levels);

% precalculate image sizes
pyramid_levels = levels+1;
width_Pyramid{1} = N;
height_Pyramid{1} = M;
for i = 2:pyramid_levels
  width_Pyramid{i}  = round(pyramid_factor_x*width_Pyramid{i-1});
  height_Pyramid{i} = round(pyramid_factor_y*height_Pyramid{i-1});
  if max(width_Pyramid{i}, height_Pyramid{i}) < startLev
    pyramid_levels = i;
    break;
  end
end

  % set up image pyramides
  I1_Pyramid = cell(pyramid_levels,1);
  I2_Pyramid = cell(pyramid_levels,1);
 
  % rescale images, apply struture texture decomposition if wanted
  for i = 1:pyramid_levels
      images = cat ( 3, imresize(I1, [height_Pyramid{i} width_Pyramid{i}], 'bicubic'), ...
                        imresize(I2, [height_Pyramid{i} width_Pyramid{i}], 'bicubic') );

    if p.stt >0
      if isfield(p, 'mask')
        mask = imresize( p.mask, [height_Pyramid{i} width_Pyramid{i}], 'nearest');
        images(:,:,1) = structure_texture_decompositionMask( images(:,:,1), p.stt, mask );
        images(:,:,2) = structure_texture_decompositionMask( images(:,:,2), p.stt, mask );
      else
        images(:,:,1) = structure_texture_decompositionMask( images(:,:,1), p.stt );
        images(:,:,2) = structure_texture_decompositionMask( images(:,:,2), p.stt );      
      end
	% hm what about this for glaciers
%    images = scale_image(double(images), 0, 1); % ? on/off?
    end
    
    I1_Pyramid{i} = images(:,:,1);
    I2_Pyramid{i} = images(:,:,2);
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% start coarse to fine processing
for level = pyramid_levels:-1:1;
  
  M = height_Pyramid{level};
  N = width_Pyramid{level};
  
  if level == pyramid_levels
    
    % initialization
    u = zeros(M,N);
    v = zeros(M,N);
    % helper fields, qq needed for TGV
    w   = zeros(M,N,4);
    pp  = zeros(M,N,4);
    qq  = zeros(M,N,8);

  else

    rescale_factor_u = width_Pyramid{level+1}/width_Pyramid{level};
    rescale_factor_v = height_Pyramid{level+1}/height_Pyramid{level};

    % prolongate to finer grid
    u = imresize(u,[M N], 'bilinear')/rescale_factor_u;
    v = imresize(v,[M N], 'bilinear')/rescale_factor_v;

    pp_tmp = pp;
    pp = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      pp(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear');
    end
    pp_tmp = qq;
    qq = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      qq(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear');
    end
    pp_tmp = w;
    w  = zeros(M,N, size( pp_tmp, 3 ));
    for i=1:size( pp_tmp, 3 )
      w(:,:,i) = imresize(pp_tmp(:,:,i), [M N], 'bilinear');
    end
    clear('pp_tmp');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if  isfield(p, 'mask')
    mask = imresize( p.mask, [height_Pyramid{level} width_Pyramid{level}], 'nearest');
  else
    mask = false (size(u));
  end
  
  if ~exist('constraints','var')
   [u, v, w, pp, qq] = tgv_run_flow(I1_Pyramid{level}, I2_Pyramid{level}, ...
                            u, v, w, pp, qq, p, mask );
  else
    
    level_constraints = constraints;rescale_cons=[0,0];
    rescale_cons(1) = width_Pyramid{level}/width_Pyramid{1};
    rescale_cons(2) = height_Pyramid{level}/height_Pyramid{1};

    level_constraints.p  = bsxfun( @times, level_constraints.p, rescale_cons);
    level_constraints.uv = bsxfun( @times, level_constraints.uv, rescale_cons);
    
    [u, v, w, pp, qq] = tgv_run_flow(I1_Pyramid{level}, I2_Pyramid{level}, ...
                            u, v, w, pp, qq, p, mask, level_constraints );
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% store final flow
flow_out = zeros(M,N,2);
flow_out(:,:,1) = u;
flow_out(:,:,2) = v;
end
