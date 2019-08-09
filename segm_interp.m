function tissue = segm_interp(dataCmp,interpSize,segm)
% segm_interp() integrates a number of tissue masks to one segmentation
% array including all masks. Missing values and overlapping tissue labels
% are interpolated using the information of nearest neighbour voxels.
%
% INPUT         dataCmp     array holding logical information about tissue
%                           types (dim1 x dim2 x dim3 x tissue)
%               interpSize  size of sphere that defines how many
%                           neighbouring voxels are included into
%                           interpolation for each missing label
%               segm        structure holding information about anatomy and
%                           masks for segmentation. This structure is
%                           computed in segm_6C.m. Only needed for plotting
%
%                                               by Jan-Ole Radecke 05/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dataTmp = dataCmp; dataTmp(dataTmp > 0) = 1;
% fixate bone compacta to ensure spongiosa is always enclosed by compacta
dataTmp2 = dataTmp(:,:,:,1:5);
compactaMask = repmat(dataTmp(:,:,:,6),[1,1,1,size(dataTmp2,4)]);
dataTmp2(compactaMask > 0) = 0; dataTmp(:,:,:,1:5) = dataTmp2;
% clean up: find ambigue labels in dataCmp, define them as missing
cleanmask = repmat(sum(dataTmp,4),[1,1,1,size(dataCmp,4)]);
dataCmp(cleanmask > 1) = 0;
kDim = size(dataCmp);
tissue = zeros(kDim(1:3));
for iCmp = 1:size(dataCmp,4)
    tissue(dataCmp(:,:,:,iCmp) > 0) = iCmp;
end % for iCmp
% checkSeg(segm,tissue);

if any(unique(tissue)' == [0:1:6]) == 0
    warning('Segmentation does not contain all relevant tissue types. Check parameters!');
end
% interpolate tissue types/missing voxels (except skull compartments, see above)
tissue(segm.headmask_h == 1 & tissue == 0) = nan;
kCount = 1;
countChange = zeros(1,6);
s = strel_bol(interpSize);
frameZeros = zeros(size(tissue)+interpSize*2);
frameZeros(interpSize+1:end-interpSize,interpSize+1:end-interpSize,interpSize+1:end-interpSize) = tissue;
checkSize = size(tissue);
tissue = frameZeros;

fprintf('Interpolation step:');
while any(any(any(isnan(tissue))))
    fprintf(' %i',kCount);
    storeTissue = tissue;
    for iX = 2:size(tissue,1)-1
        for iY = 2:size(tissue,2)-1
            for iZ = 2:size(tissue,3)-1
                if isnan(tissue(iX,iY,iZ))                                    
                    % find dominant neighbour ~= 0 for each missing element in dataTmp.tissue
                    tempTissue = tissue(iX-interpSize:iX+interpSize,iY-interpSize:iY+interpSize,iZ-interpSize:iZ+interpSize);
                    tempTissue(s == 0) = 0; % only include values surrounding the current voxel
                    if ~(sum(sum(sum(isnan(tempTissue) | tempTissue == 0))) == numel(tempTissue)) % only proceed if any element is ~= Nan | 0
                        uni = unique(tempTissue(tempTissue ~= 0 & ~isnan(tempTissue)));
                        kMax = zeros(size(uni));
                        for iUni = 1:size(uni,1)
                            kMax(iUni) = sum(sum(sum(ismember(tempTissue,uni(iUni)))));
                        end
                        maxUni = uni(kMax == (max(kMax)));
                        if sum(ismember(maxUni,3))
                            storeTissue(iX,iY,iZ) = 3; % CSF is favored if part of the maximum neighbour
                            countChange(3) = countChange(3) + 1;
                        else
                            storeTissue(iX,iY,iZ) = maxUni(1); % take first largest values if more than one tissues are max
                            countChange(maxUni(1)) = countChange(maxUni(1)) + 1;
                        end  
                    end % if any element ~= Nan | 0
                end % if isnan()   
            end % for iZ
        end % for iY
    end % for iX
    tissue = storeTissue;
    kCount = kCount + 1;
end % while any(NaN)
fprintf('\nFinished interpolation...');
tissue = tissue(interpSize+1:end-interpSize,interpSize+1:end-interpSize,interpSize+1:end-interpSize);
if ~sum(checkSize == size(tissue)) % sanity check of zero-pad removal 
    error('Wrong size of dataTmp.tissue matrix after resizing to original matrix!');
end


end

