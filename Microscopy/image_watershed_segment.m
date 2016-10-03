%
% Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
% Last Modified: Sept 30, 2016
% Use watershed tranform to detect individual cells 
% Tested on MATLAB R2011a
%

function [labelled_cells, labelled_borders] = image_watershed_segment(path, neiter)

    I = imread(path);
    if size(I,3) == 3
        I = rgb2gray(I);
    end
	
    % sobel derivative
	hy = fspecial('sobel');
	hx = hy';
	Iy = imfilter(double(I), hy, 'replicate');
	Ix = imfilter(double(I), hx, 'replicate');
	gradmag = sqrt(Ix.^2 + Iy.^2);
    
    % find foreground markers
    [foreground, ~] = image_morphological_segment(path);
    sedisk = strel('disk', 4);
    
    % erode 6 times for MIAPaCa_3.tif and 10 times for MIAPaCa_6.JPG
    for i = 1 : neiter
        foreground = imerode(foreground, sedisk);
    end
  
    foreground = bwareaopen(foreground, 10);
    I2 = I;
    I2(foreground) = 255;
    
    % find background markers
    D = bwdist(foreground);
    DL = watershed(D);
    bgm = DL == 0;
    I3 = I;
    I3(bgm) = 255;
    
    % watershed transform
    gradmag2 = imimposemin(gradmag, bgm | foreground);
    L = watershed(gradmag2);
    I4 = I;
    I4(imdilate(L == 0, ones(3, 3))) = 255;
    
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    
    if (usejava('desktop') == 1)
        figure
        subplot(2,2,1), imshow(I2), title('Foreground markers')
        subplot(2,2,2), imshow(I3), title('Watershed ridge lines')
        subplot(2,2,3), imshow(I4), title('Watershed boundary')
        subplot(2,2,4), imshow(Lrgb), title('Colored watershed label')
    end
    
    se = strel('disk', 15);
    
    % opening-by-reconstruction (Iobr)
    Ie = imerode(I, se);
    Iobr = imreconstruct(Ie, I);
    
    % opening-closing by reconstruction (Iobrcbr)
    Iobrd = imdilate(Iobr, se);
    Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
    Iobrcbr = imcomplement(Iobrcbr);
    
    % regional maxima of opening-closing by reconstruction (fgm)
    fgm = zeros(size(Iobrcbr,1), size(Iobrcbr,2));
    if size(imread(path),3) == 3
        % fgm = imregionalmin(Iobrcbr);
        fgm(Iobrcbr < prctile(reshape(Iobrcbr,[1,numel(Iobrcbr)]), 10)) = 1;
    else
        % fgm = imregionalmax(Iobrcbr);
        fgm(Iobrcbr > prctile(reshape(Iobrcbr,[1,numel(Iobrcbr)]), 90)) = 1;
    end
    
    % shrink and remove stray pixels
    se2 = strel(ones(5,5));
    fgm2 = imclose(fgm, se2);
    fgm3 = imerode(fgm2, se2);
    fgm4 = bwareaopen(fgm3, 1000);
    
    % clear border
    fgm4 = imclearborder(fgm4, 4);
    I3 = I;
    I3(fgm4) = 255;
    
    % calculate watershed ridge lines
    D = bwdist(fgm4);
    DL = watershed(D);
    bgm = DL == 0;
    
    gradmag2 = imimposemin(gradmag, bgm | fgm4);
    L = watershed(gradmag2);
    
    % create labelled matrices for identified cells and borders
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
    
    result = I;
    result(imdilate(L == 0, ones(3, 3))) = 255;
    
    if (usejava('desktop') == 1)
        figure
        subplot(2,3,1), imshow(Iobrcbr), title('Opening-closing by reconstruction')
        subplot(2,3,2), imshow(fgm), title('Foreground markers')
        subplot(2,3,3), imshow(fgm3), title('Foreground thinning')
        subplot(2,3,4), imshow(I3), title('Foreground on image')
        subplot(2,3,5), imshow(bgm), title('Watershed ridge lines')
        subplot(2,3,6), imshow(result), title('Colored watershed label')
    end
    
end
