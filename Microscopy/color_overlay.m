%
% Create a color composite image from 3 channels of microscope
% Author: Dhananjay Bhaskar
% Last Modified: April 4, 2016 
%

phaseImgs = dir(strcat('LUT-Gray',filesep,'*.tif'));
mCherryImgs = dir(strcat('LUT-Red',filesep,'*.tif'));
GFPImgs = dir(strcat('LUT-Green',filesep,'*.tif'));

if (length(phaseImgs) ~= length(mCherryImgs)) || (length(phaseImgs) ~= length(GFPImgs))
	disp('Number of images does not match.');
	return;
end

% For newer versions:
% v = VideoWriter('color_composite','Archival');
% open(v);
v = avifile('color_composite.avi', 'fps',30, 'quality',100);
mkdir('LUT-Composite')

for i = 1 : length(phaseImgs)

	[pathstr,name,ext] = fileparts(strcat('LUT-Gray',filesep,phaseImgs(i).name));
	Igray = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	[pathstr,name,ext] = fileparts(strcat('LUT-Green',filesep,GFPImgs(i).name));
	Igreen = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	[pathstr,name,ext] = fileparts(strcat('LUT-Red',filesep,mCherryImgs(i).name));
	Ired = imread(strcat(pathstr,filesep,name,ext),ext(end-2:end));
	
	% Enhance gray channel and rescale intensity
	Igraytmp = imadjust(Igray);
	Igraytmp = medfilt2(Igraytmp);
	graybackground = imopen(Igray,strel('disk',10));
	Igrayenhanced = Igraytmp - graybackground;
	
	% Enhance green and red channel
	Igenhanced = medfilt2(imadjust(Igreen));
	Irenhanced = medfilt2(imadjust(Ired));
	Iphasegreen = Igrayenhanced + Igenhanced;
	Iphasered = Igrayenhanced + Irenhanced;
	Iredenhanced = Iphasered - Igenhanced;
	Igreenenhanced = Iphasegreen - Irenhanced;
	
	output = cat(3, Iredenhanced, Igreenenhanced, zeros(size(Igreenenhanced)));
	output_name = regexprep(regexprep(strcat(name,ext),'_ch\d+',''),'-\d+-\d+','');
	imwrite(output,strcat('LUT-Composite',filesep,output_name));
	
	v = addframe(v,output); 	% For newer versions: writeVideo(v,output);
	 
end

v = close(v);
