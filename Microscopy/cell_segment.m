%
% Create a coloured composite image by combining phase, red and green fluorescent channels
% Enhance flourescent images for cell segmentation using CellProfiler
% Tested on MATLAB R2010a
% Author: Dhananjay Bhaskar
% Last Modified: April 6, 2016 
%

phaseImgs = dir(strcat('LUT-Gray', filesep, '*.tif'));
mCherryImgs = dir(strcat('LUT-Red', filesep, '*.tif'));
GFPImgs = dir(strcat('LUT-Green', filesep, '*.tif'));

if (length(phaseImgs) ~= length(mCherryImgs)) || (length(phaseImgs) ~= length(GFPImgs))
	disp('Number of images does not match.');
	return;
end

% For newer versions:
% v = VideoWriter('color_composite','Archival');
% open(v);

v = avifile('color_composite.avi', 'fps',30, 'quality',100);
v1 = avifile('LUT_Red_Enhanced.avi', 'fps', 30, 'quality', 100);
v2 = avifile('LUT_Green_Enhanced.avi', 'fps', 30, 'quality', 100);

mkdir('LUT-Composite')
mkdir('LUT-Red-Enhanced')
mkdir('LUT-Green-Enhanced')

for i = 1 : length(phaseImgs)

	[pathstr,name,ext] = fileparts(strcat('LUT-Gray',filesep,phaseImgs(i).name));
	Igray = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	[pathstr,name,ext] = fileparts(strcat('LUT-Green',filesep,GFPImgs(i).name));
	Igreen = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	[pathstr,name,ext] = fileparts(strcat('LUT-Red',filesep,mCherryImgs(i).name));
	Ired = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	
	% Make cell boundaries apparent
	Igraytmp = imadjust(Igray);
	Igraytmp = medfilt2(Igraytmp);
	graybackground = imopen(Igray, strel('disk',10));
	Igrayboundary = Igraytmp - graybackground;
	
	% Reduce noise in green and red channels
	Igenhanced = medfilt2(imadjust(Igreen));
	Irenhanced = medfilt2(imadjust(Ired));
	
	% Create composite
	Iphasegreen = Igrayboundary + Igenhanced;
	Iphasered = Igrayboundary + Irenhanced;
	Iredenhanced = Iphasered - Igenhanced;
	Igreenenhanced = Iphasegreen - Irenhanced;
	output = cat(3, Iredenhanced, Igreenenhanced, zeros(size(Igreenenhanced)));
	
	% Re-scale intensity
	rescaleintensity = 1.5;
	Igrayboundary = 2*Igrayboundary;
	
	% Enhance images for segmentation
	Iredwithboundary = Irenhanced + Igrayboundary;
	Iseggreen = rescaleintensity.*Igenhanced - Iredwithboundary;
	Igreenwithboundary = Igenhanced + Igrayboundary;
	Isegred = rescaleintensity.*Irenhanced - Igreenwithboundary;
	
	output_name = regexprep(regexprep(strcat(name, ext), '_ch\d+', ''), '-\d+-\d+','');
	
	imwrite(output, strcat('LUT-Composite', filesep, output_name));
	imwrite(Isegred, strcat('LUT-Red-Enhanced', filesep, Igreenfname, '.png'));
	imwrite(Iseggreen, strcat('LUT-Green-Enhanced', filesep, Iredfname, '.png'));
	
	% For newer versions: 
	% writeVideo(v, output);
	
	v = addframe(v, output);
	v1 = addframe(v1, Isegred);
	v2 = addframe(v2, Iseggreen);
	
end

v = close(v);
v1 = close(v1);
v2 = close(v2);
