%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% This file is a part of the MRIViewer project:                         %
% An interactive user interface for viewing and analyzing 3D structural %
% magnetic resonance images                                             %
% Copyright (C) 2016, Nima Alamatsaz, Elijah Shifnadel                  %
% All rights reserved.                                                  %
% Email: nialamat@gmail.com, es224@njit.edu                             %
% Web:   http://github.com/nalamat/mriviewer                            %
%                                                                       %
% MRIViewer is free software: you can redistribute it and/or modify     %
% it under the terms of the GNU General Public License as published by  %
% the Free Software Foundation, either version 3 of the License, or     %
% any later version.                                                    %
%                                                                       %
% MRIViewer is distributed in the hope that it will be useful,          %
% but WITHOUT ANY WARRANTY; without even the implied warranty of        %
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          %
% GNU General Public License for more details.                          %
%                                                                       %
% You should have received a copy of the GNU General Public License     %
% along with MRIViwer. If not, see <http://www.gnu.org/licenses/>.      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = imadjust3(img)
	%Parse inputs and initialize variables
	[img,imageType,lowIn,highIn,lowOut,highOut,gamma] = parseInputs(img);

	validateLowHigh(lowIn,highIn,lowOut,highOut);
	gamma = validateGamma(gamma,imageType);

	out = adjustGrayscaleImage(img,lowIn,highIn,lowOut,highOut,gamma);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = adjustGrayscaleImage(img,lIn,hIn,lOut,hOut,g)
	expansionFactor = 1;
	out = adjustArray(img, lIn, hIn, lOut, hOut, g, expansionFactor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = adjustArray(img,lIn,hIn,lOut,hOut,g,d)

%make sure img is in the range [lIn;hIn]
img(:) =  max(lIn(d,:), min(hIn(d,:),img));

out = ( (img - lIn(d,:)) ./ (hIn(d,:) - lIn(d,:)) ) .^ (g(d,:));
out(:) = out .* (hOut(d,:) - lOut(d,:)) + lOut(d,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img,imageType,low_in,high_in,low_out,high_out,gamma] = ...
	parseInputs(img)

	% Default values
	lowhigh_in  = [0; 1];
	lowhigh_out = [0; 1];
	gamma = 1;

	% IMADJUST(I)
	if ~ismatrix(img)
		error(message('images:imadjust:oneArgOnlyGrayscale','IMADJUST(I)'))
	end
	validateattributes(img, {'double' 'uint8' 'uint16' 'int16' 'single'}, ...
		{'2d'}, mfilename, 'I', 1);

	% If a user passes in a m-by-3 double array, assume it is an intensity
	% image (there is really no way to tell).
	imageType = 'intensity';

	% Turn off warning 'images:imhistc:inputHasNans' before calling STRETCHLIM and
	% restore afterwards. STRETCHLIM calls IMHIST/IMHISTC and the warning confuses
	% a user who calls IMADJUST with NaNs.
	s = warning('off','images:imhistc:inputHasNaNs');
	lowhigh_in = stretchlim(img);
	warning(s) 


	[low_in, high_in]   = splitRange(lowhigh_in, imageType);
	[low_out, high_out] = splitRange(lowhigh_out, imageType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rangeMin, rangeMax] = splitRange(range, imageType)

if numel(range) == 2
    if strcmp(imageType, 'intensity')
        rangeMin = range(1);
        rangeMax = range(2);
    else   
        % Create triples for RGB image or Colormap
        rangeMin = range(1) * ones(1,3);
        rangeMax = range(2) * ones(1,3);
    end
else
    % range is a 2 by 3 array
    rangeMin = range(1,:);
    rangeMax = range(2,:);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function validateLowHigh(lowIn,highIn,lowOut,highOut)

if any(lowIn >= highIn)
    error(message('images:imadjust:lowMustBeSmallerThanHigh'))
end

if isInvalidRange(lowIn) || isInvalidRange(highIn) ...
        || isInvalidRange(lowOut) || isInvalidRange(highOut)
    error(message('images:imadjust:parametersAreOutOfRange'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function isInvalid = isInvalidRange(range)

isInvalid = min(range) < 0 || max(range) > 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gamma = validateGamma(gamma,image_type)

if strcmp(image_type,'intensity')
    validateattributes(gamma,{'double'},{'scalar', 'nonnegative'}, ...
        mfilename, 'GAMMA', 4)
else
    validateattributes(gamma,{'double'},{'nonnegative','2d'},...
        mfilename, 'GAMMA', 4)
    if numel(gamma) == 1,
        gamma = gamma*ones(1,3);
    elseif numel(gamma) ~=3,
        error(message('images:imadjust:invalidGamma'));
    end
end
