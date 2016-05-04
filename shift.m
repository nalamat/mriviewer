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

function y=shift(x, rank, thresh)
	x = double(x);
	minx = min(x(:)); maxx = max(x(:));
	x = (x-minx)/(maxx-minx);
	y = 1./(1+exp(-rank*(x-thresh)));
	miny = min(y(:)); maxy = max(y(:));
	y = (y-miny)/(maxy-miny);
	y = y*(maxx-minx)+minx;