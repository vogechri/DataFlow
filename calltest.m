%
% wrapper calls TGV_flow -- all parameters are described there.
%
% test example: 
% calltestTGV(191, 1.25/255, 12, 3, 0.9, 10, 2, 1, 0, 0.5, 16, 0, 1 )
% calltestTGV(110, 1.25/255, 12, 3, 0.9, 10, 2, 1, 0, 0.5, 16, 0, 1, './test' )
function calltest( inr, cEps, lambda, warps, pyramid_factor, innerIts, ...
  ring, dataTerm, doTV, stt, startResolution, edgeFilter, ...
  medFilt, sFolder, mask, constraints)

addpath('./KittiIO/');
path(path,'../devkit/matlab/');
path(path,'./mex/');
doTesting = 0;
p.imageName = sprintf('KITTI_%03d', inr );


[I1, I2, flowGT, flowGT_noc, p.imageName] = ...
...  loadKITTIImage(inr, '../../data/data_stereo_flow/training/',  10, 1);
  loadKITTIImage(inr, './data/data_stereo_flow/training/',  10, 1);

flow = Data_flow(cEps, lambda, warps, pyramid_factor, I1, I2, innerIts, ...
  ring, dataTerm, doTV, stt, startResolution, ...
  edgeFilter, medFilt );

err2f = flow_error(flowGT(:,:,1:end), flow, 2);
err3f = flow_error(flowGT(:,:,1:end), flow, 3);
err4f = flow_error(flowGT(:,:,1:end), flow, 4);
err5f = flow_error(flowGT(:,:,1:end), flow, 5);

err2fn = flow_error(flowGT_noc(:,:,1:end), flow, 2);
err3fn = flow_error(flowGT_noc(:,:,1:end), flow, 3);
err4fn = flow_error(flowGT_noc(:,:,1:end), flow, 4);
err5fn = flow_error(flowGT_noc(:,:,1:end), flow, 5);

epeErr  = getEndPointError(cat(3, flow, ones(size(I1))), flowGT);
epeErrN = getEndPointError(cat(3, flow, ones(size(I1))), flowGT_noc);


% kittiStr = sprintf('DispPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', err2f, err3f, err4f, err5f, err2f, err3f, err4f, err5f);
% kittiStr = sprintf('%s\nDispPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', kittiStr, err2fn, err3fn, err4fn, err5fn, err2fn, err3fn, err4fn, err5fn);
% kittiStr = sprintf('%s\nDispEPE %.3f & %.3f\nFlowEPE %.3f & %.3f\n', kittiStr, epeErr, epeErrN, epeErr, epeErrN);

kittiStr = sprintf('FlowPix-occ 2/3/4/5 %.3f & %.3f & %.3f & %.3f\nFlowPix-noc 2/3/4/5 %.3f & %.3f & %.3f & %.3f\n', err2f, err3f, err4f, err5f, err2fn, err3fn, err4fn, err5fn);
kittiStr = sprintf('%s\nEPE %.3f & EPE(noc) %.3f\n', kittiStr, epeErr, epeErrN);
    
fprintf(kittiStr);

%   flow_write(cat(3, flow, zeros(size(flow,1), size(flow,2))), sprintf('%s%s.png', sFolder, imageName) );
%   flow_read reads it back

if exist('sFolder','var')
  if ~exist(sFolder, 'dir')
    mkdir(sFolder);
  end
  fid = fopen(sprintf('%s/RESULTS_%03d_%s.txt', sFolder, inr, date), 'w', 'n');
  if (fid ~= -1)
    fwrite(fid, kittiStr, 'char');
    fclose(fid);
  end
end

%  fprintf(1, 'Elapsed total time: %f\n', totalTime);
