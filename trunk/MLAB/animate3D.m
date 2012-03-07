%{
AUTHOR: Michael Folkerts
EMAIL: mfolkerts@physics.ucsd.edu
WEBSITE: http://bit.ly/folkerts
DATE: Feb. 2012

README: This script takes a text file as input and plots the first 3 columns
as x, y, and z coordinates then animates consecutive frames. With a little
modification, one could use a 4th column for formatting the plot shape/color.

REVISIONS:
3/6/2012 - added Octave compatability (requires ffmpeg)
3/5/2012 - scrips is now a function

%}

function animate3D(filePath, numPoints, axisMax)
	% set below to true to enable Octave compatability (requires ffmpeg)
	Use_Octave = false; %true;

	%DEFAULT PARAMETERS (could easily add more parameters to to function)

	%the following are passed as function parameters instead
	%filePath = 'data4D.txt'; % output from Aarseth code (or anything you want)
	%numPoints = 10000; % number of points in simulation data (or in each frame)
	%axisMax = 10; % range for each plot axis = -axisMax to +axisMax (in your length units)

	height = 512; % height of plot window
	width = 512; % width of plot window
	frameRate = 10;  % frame rate (frames/seconds)
	quality = 75;    % video quality (1-100)
	aviPrefix = 'plotAnimation'; % '.avi' is automatically appended

	% ADVANCED PARAMETERS (don't edit unless you know what you are doing)
	dataFormat = '%f %f %f %f'; % specifies each line to contain 4 floats
	vectorSize = 4;    % needs be at least 3, extra colums get ignored
					   % (extra columns could be used to format the data
					   %  points if you modify the code)
	dpi = 100; % used for octave print function, sets relative size of fonts etc.
	tempDir = 'tempPNG';
	ffmpegExt = '.mp4'; % change this as needed to match codecs included in your installation

	%--START SCRIPT--

	% open the file:
	FileID = fopen(filePath,'r');

	% read the file (assume x y z glyph):
	[data, elements] = fscanf(FileID,dataFormat,[vectorSize,Inf]);

	%TODO: input error checking!

	numFrames = elements/vectorSize/numPoints; % number of frames in data

	% reshape data array (coords,point,frame):
	data = reshape(data,vectorSize,numPoints,elements/vectorSize/numPoints);

	if Use_Octave
		% create temporary folder for .png files:
		mkdir(tempDir);
	else
		% initialize video object (MATLAB ONLY):
		myVideo = VideoWriter(aviPrefix); % 'Motion JPEG AVI' is default
		myVideo.FrameRate = frameRate;
		myVideo.Quality = quality;
		% open video object:
		open(myVideo);
	end

	figHandle = figure;
	figPos = get(figHandle,'position');
	figPos(3) = width;
	figPos(4) = height;
	set(figHandle,'position',figPos); 

	%the above gets ignored when printing figures in Octave (need to use printer properties instead):
	set(gcf, 'PaperPosition', [0 0 (width/dpi) (height/dpi)]);

	for frame = 1:numFrames

		if Use_Octave
			plot3(data(1,:,frame), data(2,:,frame), data(3,:,frame),'o') % scatter3 is way to slow in current version of Octave
		else
			scatter3(data(1,:,frame), data(2,:,frame), data(3,:,frame),...
				20, 'o',... % markersize and type
				'lineWidth',1, 'markerFaceColor',[.8,.8,.8] ) %grey


		end
		%set axis limits
		xlim([-axisMax axisMax])
		ylim([-axisMax axisMax])
		zlim([-axisMax axisMax])

		%set axis labels
		xlabel('x')
		ylabel('y')
		zlabel('z')
		
		pause(.2) % may not be needed, but just to be safe

		if Use_Octave
			% save plot frame as png
			print(gcf,'-dpng',['-r' int2str(dpi)],[tempDir '/plot_' int2str(frame) '.png']);
		else
			% write frame to video(MATLAB ONLY):
			writeVideo(myVideo, getframe(figHandle));
		end

	end

	if Use_Octave
		% need ffmpeg installed
		system(['ffmpeg -f image2 -qcomp ' num2str(quality/100.0) ' -r ' int2str(frameRate) ' -i ' tempDir '/plot_%d.png ' aviPrefix ffmpegExt], '-echo')
		fprintf('-------\nDone! Video has been saved as %s%s in %s \n-------\n',aviPrefix,ffmpegExt,pwd)
	else
		% close video object and finalize file(MATLAB ONLY):
		close(myVideo);
		msgbox(['Done! Video has been saved as ' myVideo.Filename ' in ' pwd])
	end

	

end
