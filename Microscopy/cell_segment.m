%
% Enhance and segment flourescent protein tagged cells for identification using CellProfiler
% Tested on MATLAB R2010a
% Author: Dhananjay Bhaskar
% Last Modified: April 6, 2016 
%

phaseImgs = dir(strcat('LUT-Gray',filesep,'*.tif'));
mCherryImgs = dir(strcat('LUT-Red',filesep,'*.tif'));
GFPImgs = dir(strcat('LUT-Green',filesep,'*.tif'));

if (length(phaseImgs) ~= length(mCherryImgs)) || (length(phaseImgs) ~= length(GFPImgs))
	disp('Number of images does not match.');
	return;
end

v1 = avifile('LUT_Red_Enhanced.avi', 'fps',30, 'quality',100);
v2 = avifile('LUT_Green_Enhanced.avi', 'fps',30, 'quality',100);
mkdir('LUT-Red-Enhanced')
mkdir('LUT-Green-Enhanced')

for i = 1 : length(phaseImgs)

	% Load images

	[pathstr,name,ext] = fileparts(strcat('LUT-Gray',filesep,phaseImgs(i).name));
	Igray = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	[pathstr,name,ext] = fileparts(strcat('LUT-Green',filesep,GFPImgs(i).name));
	Igreenfname = name;
	Igreen = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	[pathstr,name,ext] = fileparts(strcat('LUT-Red',filesep,mCherryImgs(i).name));
	Iredfname = name;
	Ired = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	

	% Identify cell boundaries
	Igraytmp = imadjust(Igray);
	Igraytmp = medfilt2(Igraytmp);
	graybackground = imopen(Igray,strel('disk',2));
	Igrayboundary = Igraytmp - graybackground;

	% Enhance green channel
	Igenhanced = medfilt2(imadjust(Igreen));

	% Enhance red channel
	Irenhanced = medfilt2(imadjust(Ired));
	
	% Re-scale intensity
	rescaleintensity = 1.5;
	Igrayboundary = 2*Igrayboundary;

	% Segment image using cell boundaries from phase channel
	Iredwithboundary = Irenhanced + Igrayboundary;
	Iseggreen = rescaleintensity.*Igenhanced - Iredwithboundary;
	Igreenwithboundary = Igenhanced + Igrayboundary;
	Isegred = rescaleintensity.*Irenhanced - Igreenwithboundary;

	output_name = regexprep(regexprep(strcat(name,ext),'_ch\d+',''),'-\d+-\d+','');
	
	imwrite(Isegred,strcat('LUT-Red-Enhanced',filesep,Iredfname,'.png'));
	imwrite(Iseggreen,strcat('LUT-Green-Enhanced',filesep,Igreenfname,'.png'));
	
	v1 = addframe(v1,Isegred);
	v2 = addframe(v2,Iseggreen);
	
end

v1 = close(v1);
v2 = close(v2);
