%
% Author: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
% Last Modified: Sept 30, 2016
% Determine outline of cancer cell lines using edge detection (derivative filter) and mathematical morphology
% Tested on MATLAB R2011a
%

function [BWeroded, BWoutline] = image_morphological_segment(path)

	I = imread(path);
    
    if size(I,3) == 3
        I = rgb2gray(I);
    end
    
	hy = fspecial('sobel');
	hx = hy';
	Iy = imfilter(double(I), hy, 'replicate');
	Ix = imfilter(double(I), hx, 'replicate');
	gradmag = sqrt(Ix.^2 + Iy.^2);
    
	[~, threshold] = edge(I, 'sobel');
	fudgeFactor = 1;
	BWs = edge(I,'sobel', threshold * fudgeFactor);
    
	se90 = strel('line', 6, 90);
	se0 = strel('line', 6, 0);
    se45 = strel('line', 6, 45);
    se135 = strel('line', 6, 135);
    
	BWsdil = imdilate(BWs, [se90 se0 se45 se135]);
	BWdfill = imfill(BWsdil, 'holes');

    % remove border objects and artifacts
	BWnobord = imclearborder(BWdfill, 4);
    BWnobord = bwareaopen(BWnobord, 1000);
    
    % smoothing
	seD = strel('diamond',1);
	BWeroded = imerode(BWnobord,seD);
	BWeroded = imerode(BWeroded,seD);

	BWoutline = bwperim(BWeroded);
	
    sedisk = strel('disk', 2);
    BWoutlinedilated = imdilate(BWoutline, sedisk);
    Segout = I;
	Segout(BWoutlinedilated) = 255;
    
    if (usejava('desktop') == 1)
        figure
        subplot(2,3,1), imshow(gradmag,[]), title('Gradient magnitude')
        subplot(2,3,2), imshow(BWs), title('Sobel filter thresholding')
        subplot(2,3,3), imshow(BWsdil), title('Dilated and filled')
        subplot(2,3,4), imshow(BWnobord), title('Clear border and small objects')
        subplot(2,3,5), imshow(BWeroded), title('Eroded binary image')
        subplot(2,3,6), imshow(Segout), title('Super-impose outline')
    end
end
