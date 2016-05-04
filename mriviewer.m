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

function varargout = mriviewer(varargin)
	% Begin initialization code - DO NOT EDIT
	gui_Singleton = 1;
	gui_State = struct('gui_Name',       mfilename, ...
					   'gui_Singleton',  gui_Singleton, ...
					   'gui_OpeningFcn', @fig_OpeningFcn, ...
					   'gui_OutputFcn',  @fig_OutputFcn, ...
					   'gui_LayoutFcn',  [] , ...
					   'gui_Callback',   []);
	if nargin && ischar(varargin{1})
		gui_State.gui_Callback = str2func(varargin{1});
	end

	if nargout
		[varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
	else
		gui_mainfcn(gui_State, varargin{:});
	end
	% End initialization code - DO NOT EDIT


% --- Outputs from this function are returned to the command line.
function varargout = fig_OutputFcn(obj, event, h)
	if isstruct(h)
		varargout{1} = h.output;
	else
		varargout{1} = -1;
	end


% --- Executes just before mriviewer is made visible.
function fig_OpeningFcn(obj, event, h, varargin)
	% Choose default command line output for mriviewer
	h.output = obj;

	% This sets up the initial plot - only do when we are invisible
	% so window can get raised using mriviewer.
	if strcmp(obj.Visible,'off')
		addlistener(h.slrSmooth, 'Value', 'PostSet', @(~,~)process_Callback(h.slrSmooth,event,h));
		addlistener(h.slrSigmoidRank, 'Value', 'PostSet', @(~,~)process_Callback(h.slrSigmoidRank,event,h));
		addlistener(h.slrSigmoidThresh, 'Value', 'PostSet', @(~,~)process_Callback(h.slrSigmoidThresh,event,h));
		
		folder = 'data';
		if (~isdir(folder))
			folder = uigetdir('', 'Select Data Folder');
			if folder == 0; delete(obj); return; end
		end
		load_folder(folder, h);
	else
		process_img(h);
	end
	
	% Update handles structure
	guidata(obj, h);


% --- Executes when user attempts to close fig.
function fig_CloseRequestFcn(obj, event, h)
	clearvars -global img_raw img_avg img sizeZ sizeY sizeX ...
		sliceZ sliceY sliceX view3D path times;
	delete(obj);


function fig_KeyPressFcn(obj, event, h)
	global view3D;
	global sizeZ sizeY sizeX;
	global sliceZ sliceY sliceX;
	
	a = gca;
	c = obj.CurrentCharacter;
	az=view3D(1); el=view3D(2);
	z=sliceZ; y=sliceY; x=sliceX;
	
	if c==28
		if     a==h.axesZ ; y=y-1;
		elseif a==h.axesY ; x=x+1;
		elseif a==h.axesX ; y=y-1;
		elseif a==h.axes3D; az=az-1; end
	elseif c==29
		if     a==h.axesZ ; y=y+1;
		elseif a==h.axesY ; x=x-1;
		elseif a==h.axesX ; y=y+1;
		elseif a==h.axes3D; az=az+1; end
	elseif c==30
		if     a==h.axesZ ; x=x-1;
		elseif a==h.axesY ; z=z+1;
		elseif a==h.axesX ; z=z+1;
		elseif a==h.axes3D; el=el+1; end
	elseif c==31
		if     a==h.axesZ ; x=x+1;
		elseif a==h.axesY ; z=z-1;
		elseif a==h.axesX ; z=z-1;
		elseif a==h.axes3D; el=el-1; end
	end
    
	z(z<1)=1; z(z>sizeZ)=sizeZ; sliceZ=z; 
	y(y<1)=1; y(y>sizeY)=sizeY; sliceY=y;
	x(x<1)=1; x(x>sizeX)=sizeX; sliceX=x;
	el(el<-90)=-90; el(el>90)=90;
	
	if a==h.axesZ || a==h.axesY || a==h.axesX
		update_views(h);
	elseif a==h.axes3D
		view3D = [az el];
		view(view3D);
	end


function btnLoad_Callback(obj, event, h)
	folder = uigetdir('', 'Select Data Folder');
	if folder ~= 0; load_folder(folder, h); end


function btnAverage_Callback(obj, event, h)
	
	process_img(h);

% --- Executes on selection change in lstFiles.
function lstFiles_Callback(obj, event, h)
	global path;
	load_file([path '/' obj.String{obj.Value}], h);


% --- Executes on button press in chkThreshold.
function process_Callback(obj, event, h)
	process_img(h);

% --- Executes on mouse press over axes background.
function axesXYZ_drag(obj, h)
	global sizeZ sizeY sizeX;
	global sliceZ sliceY sliceX;
	
	if strcmp(h.fig.SelectionType, 'normal')==0; return; end
	
    if obj == h.axesZ
        z = sliceZ;
        y = obj.CurrentPoint(1,1);
        x = obj.CurrentPoint(1,2);
    elseif obj == h.axesY
        x = -obj.CurrentPoint(1,1)+sizeX+1;
        y = sliceY;
        z = -obj.CurrentPoint(1,2)+sizeZ+1;
    elseif obj == h.axesX
        x = sliceX;
        y = obj.CurrentPoint(1,1);
        z = -obj.CurrentPoint(1,2)+sizeZ+1;
    end
    
	z(z<1)=1; z(z>sizeZ)=sizeZ; sliceZ=z; 
	y(y<1)=1; y(y>sizeY)=sizeY; sliceY=y;
	x(x<1)=1; x(x>sizeX)=sizeX; sliceX=x;
	
	update_views(h);


% --- Executes on mouse press over axes background.
function axes3D_drag(obj, h)
	global view3D;
	view3D = obj.UserData.v - (h.fig.CurrentPoint - obj.UserData.p);
	if view3D(2)<-90; view3D(2)=-90; end
	if view3D(2)> 90; view3D(2)= 90; end
	view(view3D);


% --- Executes on button press in chkSliceX, chkSliceY and chkSliceZ.
function chkSlices_Callback(obj, event, h)
	update_views(h);


function load_folder(folder, h)
	global path;
	
	files = dir([folder '/*.nii']);
	files = {files(~[files(:).isdir]).name};
	
	if isempty(files)
		msgbox('Selected folder does not contain any NIfTI (.nii) files');
		return;
	end
	
	path = folder;
	h.lstFiles.Value = 1;
	h.lstFiles.String = files;

	process_avg(h);
	load_file([path '/' h.lstFiles.String{h.lstFiles.Value}], h);


function load_file(file, h)
	global img_raw view3D;
	global sizeZ sizeY sizeX;
	global sliceZ sliceY sliceX;
	
	try
		nii = load_nii(file);
		img_raw = double(nii.img);
		% Shift and scale image data to [0 1]
		img_raw = img_raw  - min(img_raw(:));
		img_raw = img_raw ./ max(img_raw(:));
	catch
		if exist('load_nii','file')==2; disp('Error loading file');
		else disp('NIfTI toolbox is not on the path'); end
		
		img_raw = zeros(101,101,101);
		img_raw(50,50,50) = 1;
	end
	
	sizeZ=size(img_raw,3); sliceZ=(sizeZ-1)/2;
	sizeY=size(img_raw,2); sliceY=(sizeY-1)/2;
	sizeX=size(img_raw,1); sliceX=(sizeX-1)/2;
	view3D = [30 30];
	
	process_img(h);


function process_avg(h)
	global img_avg path;
	
	img_avg = zeros(91, 109, 91);
	n = 0;
	
	for i=1:length(h.lstFiles.String)
		try nii = load_nii([path '/' h.lstFiles.String{i}]);
		catch continue; end
		
		im = double(nii.img);
		if size(im)~=size(img_avg); continue; end
		
		im=im-min(im(:)); im=im./max(im(:));
		img_avg = img_avg + im;
		n = n + 1;
	end
	
	img_avg = img_avg / n;


function process_img(h)
	global img_raw img img_avg;
	global sizeZ sizeY sizeX;
	global sliceZ sliceY sliceX;
	
	set(gcf, 'pointer', 'watch');
	drawnow;
	
	set(gcf, 'pointer', 'arrow')
	% Show averaged image
	if h.chkAverage.Value == 1
		img = img_avg;
	else
		img = img_raw;
	end
	
	% Contrast enhancement
	if h.chkContrast.Value == 1
		% TODO: implement 3D contrast enhancement
		for z=1:sizeZ
			img(:,:,z) = imadjust(img(:,:,z));
		end
	end
	
	% Smooth
	h.txtSmooth.String = ['Smooth: ' num2str(round(h.slrSmooth.Value,1))];
	if h.slrSmooth.Value >= 0.1
		img = imgaussfilt3(img, h.slrSmooth.Value);
	end
	
	% Sigmoid function
	h.txtSigmoidRank.String = ['Sigmoid Rank: ' num2str(floor(h.slrSigmoidRank.Value))];
	h.txtSigmoidThresh.String = ['Sigmoid Thresh: ' num2str(round(h.slrSigmoidThresh.Value,2))];
	if h.slrSigmoidRank.Value >= 1
		img = shift(img, h.slrSigmoidRank.Value, h.slrSigmoidThresh.Value);
	end
	
	h.txtSigma1.String = ['Negative Thresh: ' num2str(round(h.slrSigma1.Value,1))];
	h.txtSigma2.String = ['Positive Thresh: ' num2str(round(h.slrSigma2.Value,1))];
	
	while true
		% Do nothing, only display raw image
		if h.rdoProcess1.Value == 1; break; end
		
		% Edge detection
		if h.rdoEdge.Value == 1
			method = h.lstEdge.String{h.lstEdge.Value};
			if strcmp(method,'Canny 3D')
				if exist('canny','file')==2; img = canny(double(img),2);
				else disp('Canny 3D toolbox is not on the path'); end
			else
				for z=1:sizeZ
					img(:,:,z) = edge(img(:,:,z),method);
				end
			end
			break;
		end
		
		% Histogram analysis for background/skull thresholding
		[a b] = histcounts(img, 40);
		[c d] = findpeaks(-a);
		thresh = b(d(1));
		[c, d] = findpeaks(a);
		d = d(b(d)>thresh);
		peak = b(d(1));
		
		% Skull masking and removal
		mask = img;
		for z=1:sizeZ
			s = mask(:,:,z);
			s = im2bw(s, thresh);
			s = bwareaopen(s, 50);
			s = imclose(s,strel('disk',3));
			s = imfill(s,'holes');
			if any(s(:))
				w = warning ('off','all');
				s = bwareafilt(s, 1);
				warning(w);
			end
			s = bwmorph(s, 'shrink', 7);
			mask(:,:,z) = s;
		end
		
		img(img<thresh) = peak;
		img = img .* mask;
		
		if h.rdoProcess2.Value == 1; break; end
		
		% Thresholding
		avg = mean(img(img~=0));
		sigma = std(img(img~=0));
		img2 = img;
		
		neglow = avg - sigma*h.slrSigma1.Value;
		neghigh = avg;
		pos = avg + sigma*h.slrSigma2.Value;
		img2(img<pos) = 0;
		%img2(neglow<=img & img<=neghigh) = 1;
		img2(pos<=img) = 1;
		img = img2;
		
		if h.rdoProcess3.Value == 1; break; end
		
		% Symmetric diffrentiation
		mask = imdilate(img,strel('disk',2));
		img = img - flip(mask,1);
		img(img<0) = 0;
		
		if h.rdoProcess4.Value == 1; break; end
		
		% Segmentation
		for z=1:sizeZ
			s = img(:,:,z);
			s = imdilate(s, strel('disk',1));
			s = imclose(s, strel('disk',2));
			s = imfill(s, 'holes');
			s = im2bw(s, thresh);
			if any(s(:))
				w = warning ('off','all');
				s = bwareafilt(s, 1);
				warning(w);
			end
			img(:,:,z) = s;
		end
		
		if h.rdoProcess5.Value == 1; break; end
		
		c = zeros(sizeZ,2);
		for z=1:sizeZ
			reg = regionprops(img(:,:,z));
			if ~isempty(reg)
				c(z,:) = reg.Centroid;
			end
		end
		
		for z=1:sizeZ
			count = 0;
			for i = -5:5
				if 1<=z+i && z+i<=sizeZ && norm(c(z,:)-c(z+i,:))<5
					count = count + 1;
				else
					count = count -1;
				end
			end
			if count < 5
				img(:,:,z) = 0;
			end
		end
		
		break;
	end
	
	update_plots(h);
	update_views(h);
	
	set(gcf, 'pointer', 'arrow')

function update_plots(h)
	global img_raw img view3D;
	global sizeZ sizeY sizeX;
	global sliceZ sliceY sliceX;
	
	axesSelected = gca;
	
	textProp = {'color','g','verticalalignment','top','horizontalalignment','right'};
	boxProp = {'xcolor',[.5 .5 .5], 'ycolor',[.5 .5 .5], 'zcolor',[.5 .5 .5]};
	
	axes(h.axes1);
	histogram(img(img~=0), 'facecolor','w','edgecolor','w', 'numbins',30);
% 	[a, b] = histcounts(img, 30);
% 	[c, d] = findpeaks(-a);
% 	minx = b(d(1)); miny = a(d(1));
% 	[c, d] = findpeaks(a);
% 	c = c(b(d)>minx); d = d(b(d)>minx);
% 	maxx = b(d(1)); maxy = a(d(1));
% 	
% 	[a, b] = histcounts(img(img>minx), 30);
% 	plot(b(1:end-1),a,'w')
	%plot(b(1:end-1),a,'w', minx, miny,'g*', maxx, maxy,'b*');
	%set(gca,'yscale','log')
	h.axes1.Color = 'k';
	h.axes1.YTick=[]; h.axes1.XTick=[]; set(h.axes1, boxProp{:});
	xl=xlim; yl=ylim; x=xl(2)-.015*diff(xl); y=yl(2)-.015*diff(yl);
	text(x,y, {'Processed Image Histogram'}, textProp{:});
	
	axes(h.axes2);
	x = 0:.01:1;
	rank=h.slrSigmoidRank.Value; rank(rank<1)=1;
	thresh=h.slrSigmoidThresh.Value; thresh(rank<1)=.5;
	plot(x,shift(x,rank,thresh),'w');
	h.axes2.Color = 'k';
	h.axes2.YTick=[]; h.axes2.XTick=[]; set(h.axes2, boxProp{:});
	xl=xlim; yl=ylim; x=xl(2)-.015*diff(xl); y=yl(2)-.015*diff(yl);
	text(x,y, {'Sigmoid Function'}, textProp{:});
	
	axes(axesSelected);

function update_views(h)
	global img img_raw view3D;
	global sizeZ sizeY sizeX;
	global sliceZ sliceY sliceX;
	
	axesSelected = gca;
	
	sliceXf = -sliceX+sizeX+1; % x flipped
	lineZ=[sliceZ sliceZ]; lineY=[sliceY sliceY]; lineX=[sliceX sliceX];
	lineZf=-lineZ+sizeZ+1; lineXf=-lineX+sizeX+1; % z and x flipped
	lineProp = {'color','g','hittest','off'};
	textProp = {'color','g','verticalalignment','top','horizontalalignment','left'};
	boxProp = {'xcolor',[.5 .5 .5], 'ycolor',[.5 .5 .5], 'zcolor',[.5 .5 .5]};
	
	axes(h.axesZ);
	sliceZr = round(sliceZ);
	imagesc(squeeze(img(:,:,sliceZr)), 'hittest','off');
	h.axesZ.YTick=[]; h.axesZ.XTick=[]; set(h.axesZ, boxProp{:});
	colormap('gray'); caxis([0 1]);
	text(2,2, {'Axial Plane' 'Top View' ['Slice ' num2str(sliceZr)]}, textProp{:});
	line(lineY,ylim,lineProp{:}); line(xlim,lineX,lineProp{:});

	axes(h.axesY);
	sliceYr = round(sliceY);
	imagesc(flip(imrotate(squeeze(img(:,sliceYr,:)),90),2), 'hittest','off');
	h.axesY.YTick=[]; h.axesY.XTick=[]; set(h.axesY, boxProp{:});
	colormap('gray'); caxis([0 1]);
	text(2,2, {'Coronal Plane' 'Front View' ['Slice ' num2str(sliceYr)]}, textProp{:});
	line(lineXf,ylim,lineProp{:}); line(xlim,lineZf,lineProp{:});

	axes(h.axesX);
	sliceXr = round(sliceX);
	imagesc(imrotate(squeeze(img(sliceXr,:,:)),90), 'hittest','off');
	h.axesX.YTick=[]; h.axesX.XTick=[]; set(h.axesX, boxProp{:});
	colormap('gray'); caxis([0 1]);
	text(2,2, {'Sagittal Plane' 'Right View' ['Slice ' num2str(sliceXr)]}, textProp{:});
	line(lineY,ylim,lineProp{:}); line(xlim,lineZf,lineProp{:});
	
	axes(h.axes3D)
	s = slice(single(flip(img,1)),sliceY(h.chkSliceY.Value==1),...
		sliceXf(h.chkSliceX.Value==1),sliceZ(h.chkSliceZ.Value==1));
	line([109 109],[1 1],[1 1]); line([1 1],[91 91],[1 1]); line([1 1],[1 1],[91 91]);
	set(s, 'edgecolor','none', 'facealpha','flat', 'hittest','off');
	for i=1:length(s); set(s(i), 'alphadata', shift(get(s(i),'cdata'),30,.08)); end
	h.axes3D.XTick=[]; h.axes3D.YTick=[]; h.axes3D.ZTick=[];
	set(h.axes3D, 'color','k', 'box','on', boxProp{:}, 'hittest','on');
	view(view3D); axis tight;
	arrow3([1 1 1], [1 1 sizeZ/4], 'color','y');
	arrow3([1 1 1], [1 sizeY/4 1], 'color','r');
	arrow3([1 1 1], [sizeX/4 1 1], 'color','b');
	
	h.axesZ.ButtonDownFcn = @(~,~)drag_start(h.axesZ,h);
	h.axesY.ButtonDownFcn = @(~,~)drag_start(h.axesY,h);
	h.axesX.ButtonDownFcn = @(~,~)drag_start(h.axesX,h);
	h.axes3D.ButtonDownFcn = @(~,~)drag_start(h.axes3D,h);
	
	axes(axesSelected);


function drag_start(obj, h)
	global view3D;
	
	if obj==h.axesZ || obj==h.axesY || obj==h.axesX
		axesXYZ_drag(obj,h);
		h.fig.WindowButtonMotionFcn = @(~,~)axesXYZ_drag(obj,h);
    elseif obj==h.axes3D
		obj.UserData = struct();
		obj.UserData.v = view3D;
		obj.UserData.p = h.fig.CurrentPoint;
		h.fig.WindowButtonMotionFcn = @(~,~)axes3D_drag(obj,h);
	end

	h.fig.WindowButtonUpFcn = @(~,~)drag_stop(obj,h);


function drag_stop(obj, h)
	h.fig.WindowButtonMotionFcn = '';
	h.fig.WindowButtonUpFcn = '';
