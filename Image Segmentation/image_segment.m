%{ 
%   Authors: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca> and MoHan Zhang <mohan_z@hotmail.com>
%   Last modified: Nov 21, 2016
%   Description: Segment cancer cell line images (Phase contrast microscopy)
%   Tested on MATLAB R2011a
%}

%% Main
function [] = image_segment()
    
	nfigs = 1;

	test_32(nfigs, 1);

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
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_6', nfigs);
		% save correct segmentations
 		save_segmentation('MIAPaCa_6_segmented.mat', [2 4 6 12 15 21 23 25 35 38 39 40 42], w_cells);
 		save_segmentation('MIAPaCa_6_segmented.mat', [16 19 20], wobrcbr_cells);
 		save_PIF('MIAPaCa_6_segmented.mat', 'MIAPaCa_6.pif');
	end
    
	function [nfigs] = test_3(nfigs, display)
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
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_3', nfigs);
		% save correct segmentations
 		save_segmentation('MIAPaCa_3_segmented.mat', [10 12 17 26], mm_fg);
 		save_segmentation('MIAPaCa_3_segmented.mat', [20 33], w_cells);
 		save_PIF('MIAPaCa_3_segmented.mat', 'MIAPaCa_3.pif');
	end

	function [nfigs] = test_32(nfigs, display)
		MIAPaCa_32 = strcat('dataset',filesep,'MIAPaCa_32.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_32, 0.9, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 9
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_32, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_32, 15, 20, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_32', nfigs);
		% save correct segmentations
		save_segmentation('MIAPaCa_32_segmented.mat', [4 6 7 12 11 13 15 17 18 21], mm_fg);
		save_segmentation('MIAPaCa_32_segmented.mat', [13 18 21 22 24], wobrcbr_cells);
		save_PIF('MIAPaCa_32_segmented.mat', 'MIAPaCa_32.pif');
    	end

	function [nfigs] = test_76(nfigs, display)
		MIAPaCa_76 = strcat('dataset',filesep,'MIAPaCa_76.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_76, 0.9, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 9
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_76, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_76, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_76', nfigs);
		% save correct segmentations
		save_segmentation('MIAPaCa_76_segmented.mat', [4 13], mm_fg);
		save_segmentation('MIAPaCa_76_segmented.mat', [5 8 9 10 12 14 21 22 24], w_cells);
		save_PIF('MIAPaCa_76_segmented.mat', 'MIAPaCa_76.pif');
	end

	function [nfigs] = test_77(nfigs, display)
		MIAPaCa_77 = strcat('dataset',filesep,'MIAPaCa_77.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_77, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 9
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_77, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_77, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_77', nfigs);
		% save correct segmentations
		save_segmentation('MIAPaCa_77_segmented.mat', [5 11], mm_fg);
		save_segmentation('MIAPaCa_77_segmented.mat', [3 4 15 16], w_cells);
		save_PIF('MIAPaCa_77_segmented.mat', 'MIAPaCa_77.pif');
	end

	function [nfigs] = test_78(nfigs, display)
		MIAPaCa_78 = strcat('dataset',filesep,'MIAPaCa_78.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_78, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 6
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_78, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_78, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_78', nfigs);
		% save correct segmentations
 		save_segmentation('MIAPaCa_78_segmented.mat', [8 11 17 24 28], mm_fg);
 		save_segmentation('MIAPaCa_78_segmented.mat', [5 11 23], w_cells);
 		save_PIF('MIAPaCa_78_segmented.mat', 'MIAPaCa_78.pif');
	end

	function [nfigs] = test_60(nfigs, display)
		MIAPaCa_60 = strcat('dataset',filesep,'MIAPaCa_60.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_60, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 6
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_60, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_60, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_60', nfigs);
		% save correct segmentations
 		save_segmentation('MIAPaCa_60_segmented.mat', [11 18 19 20 23 33 36], w_cells);
 		save_PIF('MIAPaCa_60_segmented.mat', 'MIAPaCa_60.pif');
	end
    
	function [nfigs] = test_61(nfigs, display)
		MIAPaCa_61 = strcat('dataset',filesep,'MIAPaCa_61.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_61, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 6
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_61, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_61, 15, 20, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_61', nfigs);
		% save correct segmentations
 		save_segmentation('MIAPaCa_61_segmented.mat', [3 6 8 24 25 29], w_cells);
 		save_PIF('MIAPaCa_61_segmented.mat', 'MIAPaCa_61.pif');
	end

	function [nfigs] = test_50(nfigs, display)
		MIAPaCa_50 = strcat('dataset',filesep,'MIAPaCa_50.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_50, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 9
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_50, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_50, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_50', nfigs);
		% save correct segmentations
        	save_segmentation('MIAPaCa_50_segmented.mat', [6 12 19], mm_fg);
        	save_segmentation('MIAPaCa_50_segmented.mat', [2 5 8 10 18], w_cells);
 		save_PIF('MIAPaCa_50_segmented.mat', 'MIAPaCa_50.pif');
	end

	function [nfigs] = test_55(nfigs, display)
		MIAPaCa_55 = strcat('dataset',filesep,'MIAPaCa_55.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_55, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 9
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_55, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_55, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_55', nfigs);
		% save correct segmentations
        	save_segmentation('MIAPaCa_55_segmented.mat', [19 26], mm_fg);
        	save_segmentation('MIAPaCa_55_segmented.mat', [2 6 8 10 15 25], w_cells);
 		save_PIF('MIAPaCa_55_segmented.mat', 'MIAPaCa_55.pif'); 
	end
	
	function [nfigs] = test_44(nfigs, display)
		MIAPaCa_44 = strcat('dataset',filesep,'MIAPaCa_44.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_44, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 9
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_44, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_44, 15, 20, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_44', nfigs);
		% save correct segmentations
        	save_segmentation('MIAPaCa_44_segmented.mat', [11 12 26 32], mm_fg);
        	save_segmentation('MIAPaCa_44_segmented.mat', [13 14 15 20 22], w_cells);
 		save_PIF('MIAPaCa_44_segmented.mat', 'MIAPaCa_44.pif');
	end

	function [nfigs] = test_48(nfigs, display)
		MIAPaCa_48 = strcat('dataset',filesep,'MIAPaCa_48.tif');
		% edge detection using mathematical morphology
		[mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_48, 1, 8, 1000, nfigs, display);
		% compute foreground markers
		foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
		foreground_markers(mm_fg > 0) = 1;
		sedisk = strel('disk', 4);
		for i = 1 : 6
			foreground_markers = imerode(foreground_markers, sedisk);
		end
		foreground_markers = bwareaopen(foreground_markers, 10);
		% marker-based watershed
		[w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_48, foreground_markers, nfigs, display);
		% opening by reconstruction and closing by reconstruction watershed
		[wobrcbr_cells, wobrcbr_borders, nfigs] = watershed_obrcbr_segment(MIAPaCa_48, 15, 10, 1000, nfigs, display);
		% plot results
		[nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, 'MIAPaCa_48', nfigs);
		% save correct segmentations
        	save_segmentation('MIAPaCa_48_segmented.mat', [12 13], mm_fg);
       		save_segmentation('MIAPaCa_48_segmented.mat', [4 5 12 14 21 22 23 30], w_cells);
 		save_PIF('MIAPaCa_48_segmented.mat', 'MIAPaCa_48.pif');     
	end
	
end

%% Plot results
function [fig_cnt] = plot_results(mm_fg, mm_outline, w_cells, w_borders, wobrcbr_cells, wobrcbr_borders, fname, fig_cnt)

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

	mm_plot = subplot(2,3,4); imshow(mm_outline_color), title('MM Outline')
	classes = numel(unique(reshape(mm_fg, [1, numel(mm_fg)])));
	for i = 1 : classes-1
		tmp = zeros(size(mm_fg,1), size(mm_fg,2));
		tmp(mm_fg == i) = 1;
		struct_array = regionprops(tmp, 'centroid');
		if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
			if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
				centroid = cat(1, struct_array.Centroid);
				txt = text(centroid(1), centroid(2), int2str(i));
				set(txt, 'fontsize', 10);
			end
		end
	end

	wm_plot = subplot(2,3,5); imshow(w_borders_color), title('Watershed (Marker) Outlines')
	classes = numel(unique(reshape(w_cells, [1, numel(w_cells)])));
	for i = 1 : classes-1
		tmp = zeros(size(w_cells,1), size(w_cells,2));
		tmp(w_cells == i) = 1;
		struct_array = regionprops(tmp, 'centroid');
		if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
			if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
				centroid = cat(1, struct_array.Centroid);
				txt = text(centroid(1), centroid(2), int2str(i));
				set(txt, 'fontsize', 10);
			end
		end
	end

	wobrcbr_plot = subplot(2,3,6); imshow(wobrcbr_borders_color), title('Watershed (OBRCBR) Outlines')
	classes = numel(unique(reshape(wobrcbr_cells, [1, numel(wobrcbr_cells)])));
	for i = 1 : classes-1
		tmp = zeros(size(wobrcbr_cells,1), size(wobrcbr_cells,2));
		tmp(wobrcbr_cells == i) = 1;
		struct_array = regionprops(tmp, 'centroid');
		if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
			if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
				centroid = cat(1, struct_array.Centroid);
				txt = text(centroid(1), centroid(2), int2str(i));
				set(txt, 'fontsize', 10);
			end
		end
	end

	% Print outlines

	fig = figure;
	set(mm_plot, 'Visible', 'off');
	mm_fig = copyobj(mm_plot, fig);
	set(mm_fig, 'Position', get(0, 'DefaultAxesPosition'));
	saveas(mm_fig, strcat(fname, '_MM', '.png'));

	fig = figure;
	set(wm_plot, 'Visible', 'off');
	wm_fig = copyobj(wm_plot, fig);
	set(wm_fig, 'Position', get(0, 'DefaultAxesPosition'));
	saveas(wm_fig, strcat(fname, '_WS', '.png'));

	fig = figure;
	set(wobrcbr_plot, 'Visible', 'off');
	wobrcbr_fig = copyobj(wobrcbr_plot, fig);
	set(wobrcbr_fig, 'Position', get(0, 'DefaultAxesPosition'));
	saveas(wobrcbr_fig, strcat(fname, '_WOBRCBR', '.png'));

	fig_cnt = fig_cnt + 1;
    
end

%% Save segmentation
function save_segmentation(filename, cell_labels, labelled_image)

	segmentations = zeros(size(labelled_image,1), size(labelled_image,2));
	max_cid = 0;

	if exist(filename, 'file')
		load(filename, 'segmentations', 'max_cid');
		assert(isequal(size(segmentations),size(labelled_image)), 'Matrix dimensions do not match.');  
	end

	num_correct_segmentations = numel(cell_labels);

	for i = 1 : num_correct_segmentations
		max_cid = max_cid + 1;
		class = cell_labels(i);
		segmentations(labelled_image == class) = max_cid;
	end

	save(filename, 'segmentations', 'max_cid');
    
end

%%  Convert MAT to PIF file
function save_PIF(filename, PIFname)

	segmentations = [];

	if exist(filename, 'file')
		load(filename, 'segmentations');
	else
		print('Error: MAT file does not exist.'); 
	end

	cell_ids = unique(nonzeros(segmentations));
	numPIFrows = nnz(segmentations);

	PIF_data = zeros(numPIFrows, 5);
	ind = 1;
	for i = 1 : length(cell_ids)
		[rows, cols] = find(segmentations == cell_ids(i));
		for cnt = 1 : length(rows)
			PIF_data(ind, 1) = cell_ids(i);
			PIF_data(ind, 2) = rows(cnt);
			PIF_data(ind, 3) = rows(cnt);
			PIF_data(ind, 4) = cols(cnt);
			PIF_data(ind, 5) = cols(cnt);
			ind = ind + 1;
		end
	end
    
    	if exist(PIFname, 'file')
		print('Error: PIF file already exists.');
    	else
		fileID = fopen(PIFname, 'w');
		cnt = 1;
		while cnt < ind
			fprintf(fileID, '%d CellU %d %d %d %d\n', PIF_data(cnt,:));
			cnt = cnt + 1;
		end
		fclose(fileID);
    	end

end
