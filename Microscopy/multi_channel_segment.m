%
% Enhance and segment fluorescent protein tagged cells
% Tested on MATLAB R2011a
% Author: Dhananjay Bhaskar
% Last Modified: Nov 14, 2016 
%

function [] = multi_channel_segment()

	phaseImgs = dir(strcat('LUT-Gray', filesep, '*.tif'));
	mCherryImgs = dir(strcat('LUT-Red', filesep, '*.tif'));
	GFPImgs = dir(strcat('LUT-Green', filesep, '*.tif'));

	if (length(phaseImgs) ~= length(mCherryImgs)) || (length(phaseImgs) ~= length(GFPImgs))
		disp('Number of images does not match.');
		return;
	end
  
	header = {'Frame', 'CellID', 'CentroidX', 'CentroidY', 'Area', 'Perimeter', ...
		      'EulerNumber', 'Extent', 'Solidity', 'EquivDiameter', ...
		      'Orientation' 'MajorAxisLength' 'MinorAxisLength' 'Eccentricity'};
    
	GFPCSVdata = [];
	mCherryCSVdata = [];

	GFPCSVfile = 'MATLAB_Green.csv';
	fid = fopen(GFPCSVfile, 'w');
	fprintf(fid, '%s,', header{1,1:end-1});
	fprintf(fid, '%s\n', header{1,end});
	fclose(fid);

	mCherryCSVfile = 'MATLAB_Red.csv';
	fid = fopen(mCherryCSVfile, 'w');
	fprintf(fid, '%s,', header{1,1:end-1});
	fprintf(fid, '%s\n', header{1,end});
	fclose(fid);

	for i = 1 : length(phaseImgs)

		% Load images
		[pathstr, name, ext] = fileparts(strcat('LUT-Gray', filesep,phaseImgs(i).name));
		Iphasename = name;
		Igray = imread(strcat(pathstr, filesep, name ,ext), ext(end-2:end));

		[pathstr, name, ext] = fileparts(strcat('LUT-Green', filesep,GFPImgs(i).name));
		Igreenfname = name;
		Igreen = imread(strcat(pathstr, filesep, name, ext), ext(end-2:end));

		[pathstr, name, ext] = fileparts(strcat('LUT-Red', filesep,mCherryImgs(i).name));
		Iredfname = name;
		Ired = imread(strcat(pathstr, filesep, name, ext), ext(end-2:end));

		% Convert to grayscale
		if size(Igray, 3) == 3
			Igray = rgb2gray(Igray);
		end

		if size(Igreen, 3) == 3
			Igreen = rgb2gray(Igreen);
		end

		if size(Ired, 3) == 3
			Ired = rgb2gray(Ired);
		end

		display = 0;
        
		[I_gray_bw, I_gray_artifacts] = binarize_phase(Igray, 15, 300, display);
		[I_red_enhanced, I_green_enhanced] = enhance_contrast(Igray, I_gray_artifacts, Igreen, Ired, Igreenfname, Iredfname, display);
		[w_seg_green, w_seg_red] = watershed_segment(I_green_enhanced, I_red_enhanced, Igreenfname, Iredfname, Iphasename, display);

		res_green = extract_features(i, w_seg_green);
		GFPCSVdata = [GFPCSVdata; res_green];
		res_red = extract_features(i, w_seg_red);
		mCherryCSVdata = [mCherryCSVdata; res_red];
		
		k_means_segment(Igray, I_gray_artifacts, Igreen, Ired, Iphasename, display);
		
		% pause before processing next image
		if display == 1
		    pause
		end

	end
    
	dlmwrite(GFPCSVfile, GFPCSVdata, '-append', 'precision', '%.3f', 'delimiter', ',');
	dlmwrite(mCherryCSVfile, mCherryCSVdata, '-append', 'precision', '%.3f', 'delimiter', ',');
    
end

%% Contrast enhancement
function [I_red_enhanced, I_green_enhanced] = enhance_contrast(Igray, I_gray_artifacts, Igreen, Ired, Igout, Irout, disp)

	if ~exist('Green-Enhanced', 'dir')
		mkdir('Green-Enhanced')
	end
	if ~exist('Red-Enhanced', 'dir')
		mkdir('Red-Enhanced')
	end
	if ~exist('Green-Eq', 'dir')
		mkdir('Green-Eq')
	end
	if ~exist('Red-Eq', 'dir')
		mkdir('Red-Eq')
	end

	Igray(I_gray_artifacts > 0) = 0;

	% Equalize
	Igeq = adapthisteq(medfilt2(Igreen));
	Ireq = adapthisteq(medfilt2(Ired));

	I_green_enhanced = Igeq - Igray;
	I_red_enhanced = Ireq - Igray;

	uniqLevels = unique(I_green_enhanced(:));
	I_green_enhanced(I_green_enhanced < prctile(uniqLevels, 10)) = 0;
	I_green_enhanced = imfill(adapthisteq(medfilt2(I_green_enhanced)));
    
	uniqLevels = unique(I_red_enhanced(:));
	I_red_enhanced(I_red_enhanced < prctile(uniqLevels, 10)) = 0;
	I_red_enhanced = imfill(adapthisteq(medfilt2(I_red_enhanced)));

	% Plot segmentation
	if (usejava('desktop') == 1 && disp == 1)
		figure
		subplot(2,2,1), imshow(Igeq), title('Green Equalized')
		subplot(2,2,2), imshow(I_green_enhanced), title('Green Enhanced')
		subplot(2,2,3), imshow(Ireq), title('Red Equalized')
		subplot(2,2,4), imshow(I_red_enhanced), title('Red Enhanced')
	end

	imwrite(I_green_enhanced, strcat('Green-Enhanced', filesep, Igout, '.png'));
	imwrite(I_red_enhanced, strcat('Red-Enhanced', filesep, Irout, '.png'));

	imwrite(Igeq, strcat('Green-Eq', filesep, Igout, '.png'));
	imwrite(Ireq, strcat('Red-Eq', filesep, Irout, '.png'));

end

%% Binarize phase
function [I_gray_bw, I_gray_artifacts] = binarize_phase(Igray, shift, minobjsize, disp)

	% Identify cell boundaries
	I_gray_bw = im2bw(Igray, graythresh(Igray));
	Igl = circshift(I_gray_bw, [0, -shift]);
	Igl(:,end-shift+1:end) = 0;
	Igr = circshift(I_gray_bw, [0, shift]);
	Igr(:,1:shift) = 0;
	Igu = circshift(I_gray_bw, [-shift, 0]);
	Igu(end-shift+1:end,:) = 0;
	Igd = circshift(I_gray_bw, [shift, 0]);
	Igd(1:shift,:) = 0;

	% Remove artifacts and rounded cells
	I_gray_circles = Igl & Igr & Igu & Igd;
	I_gray_circles = bwareaopen(I_gray_circles, minobjsize);
	I_gray_artifacts = imfill(I_gray_circles, 'holes');
	I_gray_bw = I_gray_bw & ~I_gray_artifacts;

	% Plot boundary
	if (usejava('desktop') == 1 && disp == 1)
		figure
		subplot(2,2,1), imshow(Igray), title('Phase channel')
		subplot(2,2,2), imshow(I_gray_bw), title('Phase binary mask')
		subplot(2,2,3), imshow(I_gray_artifacts), title('Phase rounded cells')
		subplot(2,2,4), imshow(I_gray_bw), title('Updated phase binary mask')
	end

end

%% Watershed segmentation
function [s_green, s_red] = watershed_segment(I_green_enhanced, I_red_enhanced, Igout, Irout, Ipout, disp)

	if ~exist('Watershed-Outline', 'dir')
		mkdir('Watershed-Outline')
	end
	if ~exist('Watershed-Green-Outline', 'dir')
		mkdir('Watershed-Green-Outline')
	end
	if ~exist('Watershed-Red-Outline', 'dir')
		mkdir('Watershed-Red-Outline')
	end

	% params
	minobjsize = 300;
	prc = 10;
	sedisk = strel('disk', 2);

	% Enhance green channel
	I_green_bw = zeros(size(I_green_enhanced,1), size(I_green_enhanced,2));
	uniqLevels = unique(I_green_enhanced(:));
	I_green_bw(I_green_enhanced > prctile(uniqLevels, prc)) = 1;
	I_green_bw = imfill(I_green_bw, 'holes');

	% Enhance red channel
	I_red_bw = zeros(size(I_red_enhanced,1), size(I_red_enhanced,2));
	uniqLevels = unique(I_red_enhanced(:));
	I_red_bw(I_red_enhanced > prctile(uniqLevels, prc)) = 1;
	I_red_bw = imfill(I_red_bw, 'holes');

	% Compute foreground
	I_green_fg = imerode(I_green_bw, sedisk);
	I_green_fg = bwareaopen(I_green_fg, minobjsize);
	I_red_fg = imerode(I_red_bw, sedisk);
	I_red_fg = bwareaopen(I_red_fg, minobjsize);
    
	% Watershed (green)
	hy = fspecial('sobel');
	hx = hy';
	Iy = imfilter(double(I_green_enhanced), hy, 'replicate');
	Ix = imfilter(double(I_green_enhanced), hx, 'replicate');
	gradmag_green = sqrt(Ix.^2 + Iy.^2);

	D = bwdist(I_green_fg);
	DL = watershed(D);
	bg_green = DL == 0;

	gradmag_green2 = imimposemin(gradmag_green, bg_green | imerode(I_green_fg, sedisk));
	L_green = watershed(gradmag_green2);

	% Watershed (red)
	Iy = imfilter(double(I_red_enhanced), hy, 'replicate');
	Ix = imfilter(double(I_red_enhanced), hx, 'replicate');
	gradmag_red = sqrt(Ix.^2 + Iy.^2);

	D = bwdist(I_red_fg);
	DL = watershed(D);
	bg_red = DL == 0;

	gradmag_red2 = imimposemin(gradmag_red, bg_red | imerode(I_red_fg, sedisk));
	L_red = watershed(gradmag_red2);

	% Plot segmentation
	if (usejava('desktop') == 1 && disp == 1)
		figure
		subplot(2,3,1), imshow(I_green_enhanced), title('Enhanced green channel')
		subplot(2,3,2), imshow(I_green_fg), title('Green foreground mask')
		subplot(2,3,3), imshow(label2rgb(L_green, 'jet', 'w', 'shuffle')), title('Green watershed')
		subplot(2,3,4), imshow(I_red_enhanced), title('Enhanced red channel')
		subplot(2,3,5), imshow(I_red_fg), title('Red foreground mask')
		subplot(2,3,6), imshow(label2rgb(L_red, 'jet', 'w', 'shuffle')), title('Red watershed')
	end
    
	% Remove largest segment (background) and compute borders
	classes = numel(unique(reshape(L_green, [1, numel(L_green)])));
	green_borders = zeros(size(L_green,1), size(L_green,2));
	L_g = zeros(size(L_green,1), size(L_green,2));
	for j = 1 : classes-1
		tmp = zeros(size(L_green,1), size(L_green,2));
		tmp(L_green == j) = 1;
		L_g(L_green == j) = j;
		green_borders(bwperim(tmp) > 0) = j;
	end
	L_g(L_g == mode(L_g(L_g~=0))) = 0;
	green_borders(green_borders == mode(green_borders(green_borders~=0))) = 0;
    
	classes = numel(unique(reshape(L_red, [1, numel(L_red)])));
	red_borders = zeros(size(L_red,1), size(L_red,2));
	L_r = zeros(size(L_red,1), size(L_red,2));
	for j = 1 : classes-1
		tmp = zeros(size(L_red,1), size(L_red,2));
		tmp(L_red == j) = 1;
		L_r(L_red == j) = j;
		red_borders(bwperim(tmp) > 0) = j;
	end
	L_r(L_r == mode(L_r(L_r~=0))) = 0;
	red_borders(red_borders == mode(red_borders(red_borders~=0))) = 0;

	% Binarize and measure region properties
	BW_green = zeros(max(size(L_g,1),size(L_r,1)), max(size(L_g,2),size(L_r,2)));
	BW_green(L_g > 0) = 1;
	BW_green = imclearborder(BW_green, 4);
	L_g(BW_green == 0) = 0;
	s_green = regionprops(L_g, 'Centroid', 'Area', 'Perimeter', 'EulerNumber', 'Extent', 'Solidity', ...
	'EquivDiameter', 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity');

	BW_red = zeros(max(size(L_g,1),size(L_r,1)), max(size(L_g,2),size(L_r,2)));
	BW_red(L_r > 0) = 1;
	BW_red = imclearborder(BW_red, 4);
	L_r(BW_red == 0) = 0;
	s_red = regionprops(L_r, 'Centroid', 'Area', 'Perimeter', 'EulerNumber', 'Extent', 'Solidity', ...
	'EquivDiameter', 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Eccentricity');
    
	% Plot overlay (ellipse fit)
	if (usejava('desktop') == 1 && disp == 1)
		figure
		imshow(BW_red | BW_green)
		hold on
		phi = linspace(0, 2*pi, 50);
		cosphi = cos(phi);
		sinphi = sin(phi);
		for k = 1 : length(s_green)
		    xbar = s_green(k).Centroid(1);
		    ybar = s_green(k).Centroid(2);
		    a = s_green(k).MajorAxisLength/2;
		    b = s_green(k).MinorAxisLength/2;
		    theta = pi*s_green(k).Orientation/180;
		    R = [ cos(theta)   sin(theta)
		         -sin(theta)   cos(theta)];
		    xy = [a*cosphi; b*sinphi];
		    xy = R*xy;
		    x = xy(1,:) + xbar;
		    y = xy(2,:) + ybar;
		    plot(x, y, 'g', 'LineWidth', 2);
		end
		for k = 1 : length(s_red)
		    xbar = s_red(k).Centroid(1);
		    ybar = s_red(k).Centroid(2);
		    a = s_red(k).MajorAxisLength/2;
		    b = s_red(k).MinorAxisLength/2;
		    theta = pi*s_red(k).Orientation/180;
		    R = [ cos(theta)   sin(theta)
		         -sin(theta)   cos(theta)];
		    xy = [a*cosphi; b*sinphi];
		    xy = R*xy;
		    x = xy(1,:) + xbar;
		    y = xy(2,:) + ybar;
		    plot(x, y, 'r', 'LineWidth', 2);
		end
		hold off
	end
    
	% Create figure with colored cell outlines (use imoverlay if available)
	hax = figure;
	red = I_green_enhanced + I_red_enhanced;
	green = I_green_enhanced + I_red_enhanced;
	blue = I_green_enhanced + I_red_enhanced;
	[I, J] = find(green_borders > 0);
	ind = sub2ind(size(I_green_enhanced), I, J);
	green(ind) = 255;
	[I, J] = find(red_borders > 0);
	ind = sub2ind(size(I_red_enhanced), I, J);
	red(ind) = 255;
	out_image = cat(3, red, green, blue);
	if (usejava('desktop') == 1 && disp == 1)
		imshow(out_image);
		print(strcat('Watershed-Outline', filesep, Ipout, '.png'), '-dpng');
	else
		set(hax, 'Visible', 'off');
		imshow(out_image);
		saveas(hax, strcat('Watershed-Outline', filesep, Ipout, '.png'));
	end

	% Create figure showing green boundaries
	if (usejava('desktop') ~= 1 || disp == 0)
		hax = figure;
		set(hax, 'Visible', 'off');
		Seg_green = I_green_enhanced;
		Seg_green(green_borders > 0) = 255;
		imshow(Seg_green);
		saveas(hax, strcat('Watershed-Green-Outline', filesep, Igout, '.png'));
	end

	% Create figure showing red boundaries
	if (usejava('desktop') ~= 1 || disp == 0)
		hax = figure;
		set(hax, 'Visible', 'off');
		Seg_red = I_red_enhanced;
		Seg_red(red_borders > 0) = 255;
		imshow(Seg_red);
		saveas(hax, strcat('Watershed-Red-Outline', filesep, Irout, '.png'));
	end
    
end

%% K-means Colour Segmentation
function [] = k_means_segment(Igray, I_gray_artifacts, Igreen, Ired, Ipout, disp)

	if ~exist('Composite', 'dir')
		mkdir('Composite')
	end

	Igray(I_gray_artifacts > 0) = 0;

	% Equalize
	Igeq = adapthisteq(medfilt2(Igreen));
	Ireq = adapthisteq(medfilt2(Ired));

	Irgb = cat(3, Ireq+Igray, Igeq+Igray, zeros(size(Igray)));

	cform = makecform('srgb2lab');
	lab_rgb = applycform(Irgb, cform);

	ab = double(lab_rgb(:,:,2:3));
	nrows = size(ab, 1);
	ncols = size(ab, 2);
	ab = reshape(ab, nrows*ncols, 2);

	nColors = 4;
	% repeat the clustering 3 times to avoid local minima
	[cluster_idx, cluster_center] = kmeans(ab, nColors, 'distance', 'sqEuclidean', 'Replicates', 3);

	pixel_labels = reshape(cluster_idx, nrows, ncols);
	segmented_images = cell(1, 3);
	rgb_label = repmat(pixel_labels, [1 1 3]);

	for k = 1 : nColors
		color = Irgb;
		color(rgb_label ~= k) = 0;
		segmented_images{k} = color;
	end

	if (usejava('desktop') == 1 && disp == 1)
		figure
		subplot(2,2,1), imshow(segmented_images{1}), title('Cluster 1')
		subplot(2,2,2), imshow(segmented_images{2}), title('Cluster 2')
		subplot(2,2,3), imshow(segmented_images{3}), title('Cluster 3')
		subplot(2,2,4), imshow(segmented_images{4}), title('Cluster 4')
	end

	imwrite(Irgb, strcat('Composite', filesep, Ipout, '.png'));

end

%% Extract features from segmentation
function [result] = extract_features(frame, seg)

	centroids = cat(1, seg.Centroid);
	area = [seg.Area];
	perim = [seg.Perimeter];
	enum = [seg.EulerNumber];
	extent = [seg.Extent];
	solidity = [seg.Solidity];
	eqDiameter = [seg.EquivDiameter];
	orient = [seg.Orientation];
	maj_axis = [seg.MajorAxisLength];
	min_axis = [seg.MinorAxisLength];
	ecc = [seg.Eccentricity];

	result = [];
	cnt = 1;

	for j = 1 : size(seg)
		if isnan(centroids(j,1))
		    continue
		end
		result(cnt,:) = [frame cnt centroids(j,1) centroids(j,2) area(j) perim(j) ...
		              enum(j) extent(j) solidity(j) eqDiameter(j) ...
		              orient(j) maj_axis(j) min_axis(j) ecc(j)];
		cnt = cnt + 1;
	end
    
end
