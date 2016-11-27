%{ 
%   Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
%   Last modified: Oct 23, 2016
%   Description: Mathematical morphology based segmentation
%   Tested on MATLAB R2011a
%}

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
