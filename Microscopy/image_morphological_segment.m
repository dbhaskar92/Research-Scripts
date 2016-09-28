%
% Determine outline of cancer cell lines using edge detection (derivative filter) and mathematical morphology 
% Author: Dhananjay Bhaskar
% Last Modified: Sept 27, 2016 
%

function [res] = image_morphological_segment(path)

	I = imread(path);
    
	I = rgb2gray(I);
	
    figure
    
	hy = fspecial('sobel');
	hx = hy';
	Iy = imfilter(double(I), hy, 'replicate');
	Ix = imfilter(double(I), hx, 'replicate');
	gradmag = sqrt(Ix.^2 + Iy.^2);
    subplot(2,3,1), imshow(gradmag,[]), title('Gradient magnitude')

	[~, threshold] = edge(I, 'sobel');
	fudgeFactor = 1;
	BWs = edge(I,'sobel', threshold * fudgeFactor);
    subplot(2,3,2), imshow(BWs), title('Sobel filter thresholding')
    
	se90 = strel('line', 6, 90);
	se0 = strel('line', 6, 0);
    se45 = strel('line', 6, 45);
    se135 = strel('line', 6, 135);

	BWsdil = imdilate(BWs, [se90 se0 se45 se135]);
	BWdfill = imfill(BWsdil, 'holes');
    subplot(2,3,3), imshow(BWsdil), title('Dilated and filled')

	BWnobord = imclearborder(BWdfill, 4);
    BWnobord = bwareaopen(BWnobord, 1000);
    subplot(2,3,4), imshow(BWnobord), title('Clear border and small objects')
    
	seD = strel('diamond',1);
	BWfinal = imerode(BWnobord,seD);
	BWfinal = imerode(BWfinal,seD);
    subplot(2,3,5), imshow(BWfinal), title('Eroded binary image')

	BWoutline = bwperim(BWfinal);
	
    sedisk = strel('disk', 2);
    BWoutlinedilated = imdilate(BWoutline, sedisk);
    Segout = I;
	Segout(BWoutlinedilated) = 255;
    subplot(2,3,6), imshow(Segout), title('Super-impose outline')
    
    res = imerode(BWfinal, sedisk);
