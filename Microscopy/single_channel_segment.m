%{ 
%   Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
%   Last modified: Oct 23, 2016
%   Description: Segment cancer cell line images (Phase contrast microscopy)
%   Tested on MATLAB R2011a
%}

%% Main
function [] = single_channel_segment()

	nfigs = 1;

	test_1(nfigs, 1);
    
	function [nfigs] = test_1(nfigs, display)
		MIAPaCa_6 = 'MIAPaCa_6.JPG';
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_6, 1, 6, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 10
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_6, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_6, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, nfigs);
		% save correct segmentations
		save_segmentation('MIAPaCa_6_segmented.mat', [2 4 6 12 15 21 23 25 35 38 39 40 42], w_cells);
		save_segmentation('MIAPaCa_6_segmented.mat', [16 19 20], wobrcbr_cells);
	end
    
	function [nfigs] = test_2(nfigs, display)
		MIAPaCa_3 = 'MIAPaCa_3.tif';
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_3, 1, 6, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 6
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_3, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_3, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, nfigs);
		% save correct segmentations
		save_segmentation('MIAPaCa_3_segmented.mat', [10 12 17 26], mm_fg);
		save_segmentation('MIAPaCa_3_segmented.mat', [20 33], w_cells);
	end
    
end

%% Mathematical Morphology Segmentation
function [labelled_cells, labelled_borders, fig_cnt] = morphological_segment(path, fudgefactor, selength, minobjsize, fig_cnt, disp)

	I = imread(path);
	if size(I, 3) == 3
		I = rgb2gray(I);
	end
    
	[~, threshold] = edge(I, 'sobel');
	BWs = edge(I, 'sobel', threshold * fudgefactor);

	se90 = strel('line', selength, 90);
	se0 = strel('line', selength, 0);
	se45 = strel('line', selength, 45);
	se135 = strel('line', selength, 135);

	% fill holes
	BWsdil = imdilate(BWs, [se90 se0 se45 se135]);
	BWdfill = imfill(BWsdil, 'holes');

	% remove border objects and artifacts
	BWnobord = imclearborder(BWdfill, 4);
	BWnobord = bwareaopen(BWnobord, minobjsize);

	% smoothing
	seD = strel('diamond', 1);
	BWeroded = imerode(BWnobord, seD);
	BWeroded = imerode(BWeroded, seD);

	BWoutline = bwperim(BWeroded);

	% dilate outline for display
	sedisk = strel('disk', 2);
	BWoutlinedilated = imdilate(BWoutline, sedisk);
	Segout = I;
	Segout(BWoutlinedilated) = 255;

	if (usejava('desktop') == 1 && disp == 1)
		figure(fig_cnt)
		subplot(2,2,1), imshow(BWs), title('Sobel Edge Thresholding')
		subplot(2,2,2), imshow(BWsdil), title('Dilated and Filled')
		subplot(2,2,3), imshow(BWnobord), title('Clear border and small objects')
		subplot(2,2,4), imshow(Segout), title('Eroded Outline')
		fig_cnt = fig_cnt + 1;
	end

	% label foreground and outline
	cc = bwconncomp(BWeroded);
	L = labelmatrix(cc);

	classes = numel(unique(reshape(L, [1, numel(L)])));
	labelled_borders = zeros(size(L,1), size(L,2));
	labelled_cells = zeros(size(L,1), size(L,2));
	for i = 1 : classes-1
		tmp = zeros(size(L,1), size(L,2));
		tmp(L == i) = 1;
		labelled_cells(L == i) = i;
		labelled_borders(bwperim(tmp) > 0) = i;
	end

end

%% Watershed segmentation using supplied foreground markers
function [labelled_cells, labelled_borders, fig_cnt] = watershed_segment(path, foreground, fig_cnt, disp)

	I = imread(path);
	if size(I,3) == 3
		I = rgb2gray(I);
	end

	% computer gradient magnitude
	hy = fspecial('sobel');
	hx = hy';
	Iy = imfilter(double(I), hy, 'replicate');
	Ix = imfilter(double(I), hx, 'replicate');
	gradmag = sqrt(Ix.^2 + Iy.^2);

	% given foreground markers
	foreground_markers = I;
	foreground_markers(foreground) = 255;

	% find background markers
	D = bwdist(foreground);
	DL = watershed(D);
	bgm = DL == 0;
	watershed_ridge = I;
	watershed_ridge(bgm) = 255;

	gradmag2 = imimposemin(gradmag, bgm | foreground);
	L = watershed(gradmag2);
	watershed_boundary = I;
	watershed_boundary(imdilate(L == 0, ones(3, 3))) = 255;

	Lrgb = label2rgb(L, 'jet', 'k', 'shuffle');

	if (usejava('desktop') == 1 && disp == 1)
		figure(fig_cnt)
		subplot(2,2,1), imshow(foreground_markers), title('Foreground markers')
		subplot(2,2,2), imshow(watershed_ridge), title('Background markers')
		subplot(2,2,3), imshow(watershed_boundary), title('Watershed boundary')
		subplot(2,2,4), imshow(Lrgb), title('Colored watershed label')
		fig_cnt = fig_cnt + 1;
	end

	classes = numel(unique(reshape(L, [1, numel(L)])));
	labelled_borders = zeros(size(L,1), size(L,2));
	labelled_cells = zeros(size(L,1), size(L,2));
	for i = 1 : classes-1
		tmp = zeros(size(L,1), size(L,2));
		tmp(L == i) = 1;
		labelled_cells(L == i) = i;
		labelled_borders(bwperim(tmp) > 0) = i;
	end
	labelled_borders(labelled_borders == mode(labelled_borders(labelled_borders~=0))) = 0;
	labelled_cells(labelled_cells == mode(labelled_cells(labelled_cells~=0))) = 0;
    
end

%% Watershed segmentation using opening and closing by reconstruction (OBRCBR)
function [labelled_cells, labelled_borders, fig_cnt] = watershed_obrcbr_segment(path, seradius, thresh_prc, minobjsize, fig_cnt, disp)

	I = imread(path);
	if size(I,3) == 3
		I = rgb2gray(I);
	end

	% computer gradient magnitude
	hy = fspecial('sobel');
	hx = hy';
	Iy = imfilter(double(I), hy, 'replicate');
	Ix = imfilter(double(I), hx, 'replicate');
	gradmag = sqrt(Ix.^2 + Iy.^2);

	% opening-by-reconstruction (Iobr)
	se = strel('disk', seradius);
	Ie = imerode(I, se);
	Iobr = imreconstruct(Ie, I);

	% opening-closing by reconstruction (Iobrcbr)
	Iobrd = imdilate(Iobr, se);
	Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
	Iobrcbr = imcomplement(Iobrcbr);

	% regional maxima of opening-closing by reconstruction (fgm)
	fgm = zeros(size(Iobrcbr,1), size(Iobrcbr,2));
	if size(imread(path),3) == 3
		%fgm = imregionalmin(Iobrcbr);
		fgm(Iobrcbr < prctile(reshape(Iobrcbr,[1,numel(Iobrcbr)]), thresh_prc)) = 1;
	else
		%fgm = imregionalmax(Iobrcbr);
		fgm(Iobrcbr > prctile(reshape(Iobrcbr,[1,numel(Iobrcbr)]), 100-thresh_prc)) = 1;
	end
    
	% shrink and remove stray pixels
	se2 = strel(ones(5,5));
	fgm2 = imclose(fgm, se2);
	fgm3 = imerode(fgm2, se2);
	fgm4 = bwareaopen(fgm3, minobjsize);

	% clear border
	fgm4 = imclearborder(fgm4, 4);

	foreground = I;
	foreground(fgm4) = 255;

	% calculate watershed ridge lines
	D = bwdist(fgm4);
	DL = watershed(D);
	bgm = DL == 0;

	gradmag2 = imimposemin(gradmag, bgm | fgm4);
	L = watershed(gradmag2);

	watershed_boundary = I;
	watershed_boundary(imdilate(L == 0, ones(3, 3))) = 255;

	Lrgb = label2rgb(L, 'jet', 'k', 'shuffle');

	if (usejava('desktop') == 1 && disp == 1)
		figure(fig_cnt)
		subplot(2,3,1), imshow(Iobrcbr), title('Opening-closing by reconstruction')
		subplot(2,3,2), imshow(fgm), title('Foreground markers')
		subplot(2,3,3), imshow(foreground), title('Foreground on image')
		subplot(2,3,4), imshow(bgm), title('Background markers')
		subplot(2,3,5), imshow(watershed_boundary), title('Watershed boundary')
		subplot(2,3,6), imshow(Lrgb), title('Colored watershed label')
		fig_cnt = fig_cnt + 1;
	end

	classes = numel(unique(reshape(L, [1, numel(L)])));
	labelled_borders = zeros(size(L,1), size(L,2));
	labelled_cells = zeros(size(L,1), size(L,2));
	for i = 1 : classes-1
		tmp = zeros(size(L,1), size(L,2));
		tmp(L == i) = 1;
		labelled_cells(L == i) = i;
		labelled_borders(bwperim(tmp) > 0) = i;
	end
	labelled_borders(labelled_borders == mode(labelled_borders(labelled_borders~=0))) = 0;
	labelled_cells(labelled_cells == mode(labelled_cells(labelled_cells~=0))) = 0;
    
end

%% Plot results
function [fig_cnt] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, fig_cnt)

	mm_fg_color = label2rgb(mm_fg, 'jet', [.7 .7 .7], 'shuffle');
	mm_outline_color = label2rgb(imdilate(mm_outline, ones(3,3)), 'jet', [.7 .7 .7], 'shuffle');
	w_cells_color = label2rgb(w_cells, 'jet', [.7 .7 .7], 'shuffle');
	w_borders_color = label2rgb(imdilate(w_borders, ones(3,3)), 'jet', [.7 .7 .7], 'shuffle');
	wobrcbr_cells_color = label2rgb(wobrcbr_cells, 'jet', [.7 .7 .7], 'shuffle');
	wobrcbr_borders_color = label2rgb(imdilate(wobrcbr_borders, ones(3,3)), 'jet', [.7 .7 .7], 'shuffle');

	figure(fig_cnt)

	subplot(2,3,1), imshow(mm_fg_color), title('MM Foreground')
	subplot(2,3,2), imshow(w_cells_color), title('Watershed (Marker) Foreground')
	subplot(2,3,3), imshow(wobrcbr_cells_color), title('Watershed (OBRCBR) Foreground')
    
	subplot(2,3,4), imshow(mm_outline_color), title('MM Outline')
	classes = numel(unique(reshape(mm_fg, [1, numel(mm_fg)])));
	for i = 1 : classes-1
		tmp = zeros(size(mm_fg,1), size(mm_fg,2));
		tmp(mm_fg == i) = 1;
		struct_array = regionprops(tmp, 'centroid');
		if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
			if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
				centroid = cat(1, struct_array.Centroid);
				txt = text(centroid(1), centroid(2), int2str(i));
				set(txt, 'fontsize', 5);
			end
		end
	end
    
	subplot(2,3,5), imshow(w_borders_color), title('Watershed (Marker) Outlines')
	classes = numel(unique(reshape(w_cells, [1, numel(w_cells)])));
	for i = 1 : classes-1
		tmp = zeros(size(w_cells,1), size(w_cells,2));
		tmp(w_cells == i) = 1;
		struct_array = regionprops(tmp, 'centroid');
		if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
			if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
				centroid = cat(1, struct_array.Centroid);
				txt = text(centroid(1), centroid(2), int2str(i));
				set(txt, 'fontsize', 5);
			end
		end
	end

	subplot(2,3,6), imshow(wobrcbr_borders_color), title('Watershed (OBRCBR) Outlines')
	classes = numel(unique(reshape(wobrcbr_cells, [1, numel(wobrcbr_cells)])));
	for i = 1 : classes-1
		tmp = zeros(size(wobrcbr_cells,1), size(wobrcbr_cells,2));
		tmp(wobrcbr_cells == i) = 1;
		struct_array = regionprops(tmp, 'centroid');
		if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
			if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
				centroid = cat(1, struct_array.Centroid);
				txt = text(centroid(1), centroid(2), int2str(i));
				set(txt, 'fontsize', 5);
			end
		end
	end

	fig_cnt = fig_cnt + 1;
    
end

%% Save segmentation
function save_segmentation(filename, cell_labels, labelled_image)

	segmentations = zeros(size(labelled_image,1), size(labelled_image,2));

	if exist(filename, 'file')
	   load(filename, 'segmentations');
	   assert(isequal(size(segmentations),size(labelled_image)), 'Matrix dimensions do not match.');  
	end

	num_correct_segmentations = numel(cell_labels);

	for i = 1 : num_correct_segmentations
		class = cell_labels(i);
		segmentations(labelled_image == class) = i;
	end

	save(filename, 'segmentations');
	
end
