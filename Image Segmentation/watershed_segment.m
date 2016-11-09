%{ 
%   Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
%   Last modified: Oct 23, 2016
%   Description: Watershed segmentation using supplied foreground markers
%   Tested on MATLAB R2011a
%}

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
