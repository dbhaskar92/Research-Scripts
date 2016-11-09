%{ 
%   Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
%   Last modified: Oct 23, 2016
%   Description: Watershed segmentation using opening and closing by reconstruction (OBRCBR)
%   Tested on MATLAB R2011a
%}

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
