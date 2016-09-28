%
% Use watershed tranform to detect individual cells 
% Author: Dhananjay Bhaskar
% Last Modified: Sept 27, 2016 
%

function image_watershed_segment(path)

    I = imread(path);
	
	I = rgb2gray(I);
	
	hy = fspecial('sobel');
	hx = hy';
	Iy = imfilter(double(I), hy, 'replicate');
	Ix = imfilter(double(I), hx, 'replicate');
	gradmag = sqrt(Ix.^2 + Iy.^2);
    
    foreground = image_morphological_segment(path);
    sedisk = strel('disk', 4);
    for i = 1:10
        foreground = imerode(foreground, sedisk);
    end
    foreground = bwareaopen(foreground, 10);
    figure
    I2 = I;
    I2(foreground) = 255;
    subplot(2,2,1), imshow(I2), title('Foreground markers')
    
    D = bwdist(foreground);
    DL = watershed(D);
    bgm = DL == 0;
    I3 = I;
    I3(bgm) = 255;
    subplot(2,2,2), imshow(I3), title('Watershed ridge lines')
    
    gradmag2 = imimposemin(gradmag, bgm | foreground);
    L = watershed(gradmag2);
    I4 = I;
    I4(imdilate(L == 0, ones(3, 3))) = 255;
	subplot(2,2,3), imshow(I4), title('Watershed boundary')
    
    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    subplot(2,2,4), imshow(Lrgb), title('Colored watershed label')
