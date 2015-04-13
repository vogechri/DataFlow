%%%   TV/TGV variational optical flow
%%%
%%% An Evaluation of Data Costs for Optical Flow, 
%%% C. Vogel, S. Roth, K. Schindler, 
%%% In: German Conference on Pattern Recognition (GCPR), Saarbrücken, Germany, 2013
%
% Content:
% optical flow computation in a pyramidal scheme with Total/
% Total generalized varation regularization of the motion fields
%
% Paramters: 
% I1,I2: input images
% cEps: TCensus data term epsilon, default: 1.25/255
% lambda: strength of the data term; higher values less regularization
% warps: the number of warping steps in the algorithm
% pyramid_factor: scale in image pyramid, default 0.9
% innerIts: the number of inner iterations performed
% ring: if data term is CSAD/Census defines th size of the patch 
%       reasonable values 1-3, dafault: 2
% dataTerm: 0:SAD, 1: Census, 2: CSAD
% doTV: 1 TV-regularization  (piecewise constant)
% doTV: 0 TGV^2-regularization (piecewise linear )
% stt : structure texture preprocessing - 0: none, 1: full (default:0.9)
%       only reasonable for values > 0.5 and with SAD data term
% startResolution: minimal image size in the pyramid, default: 16
% edgeFilter: weaker smoothness along edges: 0: off 1: on
% medFilt: post warp median Filter -- 0: off, 1: peaks, 2: all
% mask: false in areas flow should be computed - if given
% constraints: point to point constraints eg. from sparse matches
%             constrainst are supposed to be pixel by 2d motion constraints
%             telling pixel p has motion vector (u,v)
%
% Output:
% flow: a 2 dimensional motion field
%
function flow = TGV_flow(cEps, lambda, warps, pyramid_factor, ...
  I1, I2, innerIts, ring, dataTerm, doTV, stt, startResolution, ...
  edgeFilter, medFilt, mask, constraints)

path(path,'./mex/');

if ~exist('startResolution', 'var')
  startResolution = 16;
end
if ~exist('doTV', 'var')
  doTV = 0;
end
if ~exist('stt', 'var')
  stt = 0;
end

p.cEps           = cEps;    % threshold for Ternary Census
p.lambda         = lambda;  % default: 12.33 strength of data term
p.maxits         = innerIts;% default: 20 number of inner iterations
p.warps          = warps;   % number of warps, the more the better, default: 5
p.pyramid_factor = pyramid_factor; % pyramid level spacing : default 0.9
p.doTV           = doTV;    % 0: total variation, 1: 1.st grade total generalized variation
p.stt            = stt;     % 0 : off, the higher the more robust to lighting variations
p.dataTerm       = dataTerm;% 0: SAD, 1: CENSUS, 2: CSAD - latter only doMatlab =0
p.ring           = 2;       % for CSAD, CENSUS - all above 3 is taking lots of memory
p.medFilt        = 1;       % 0: off, 1: peaks, 2: all
p.doStereo       = 0;       % 1: stereo only ie. disparities assuming rectified images
p.edgeFilter     = 0;       % weaker smoothness along image edges : 0: off, 1: on

% the mask defines the area of interesst AOI. pixels marked with true 
% are not considered in the evaluation preventing motion estimation in
% the specified areas. These can e.g. be areas which are known to be static.
% The benefits are: 
% 1. faster computation
% 2. Preventing a bad data term in these areas to deliver 
%    misleading/wrong flow estiamtions.
if exist('mask','var')
  p.mask = mask;
end
if exist('ring','var')
  p.ring = ring;
end
if exist('medFilt','var')
  p.medFilt = medFilt;
end
if exist('edgeFilter','var')
  p.edgeFilter = edgeFilter;
end

% print parameters
p

avail=license('checkout','Image_Toolbox');
while avail == 0
  %     disp('License not available, waiting...');
  pause(1.0);
  avail=license('checkout','Image_Toolbox');
end

  ticID = tic;
if exist('constraints','var')
  flow = pyramidFlow(I1, I2, p, startResolution, constraints);
else
  flow = pyramidFlow(I1, I2, p, startResolution );  
end
  totalTime = toc(ticID);

  fprintf(1, 'Elapsed total time: %f\n', totalTime);
end
