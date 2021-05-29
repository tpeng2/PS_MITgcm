%convert png files to a movie
%% setup
workingDir = path_output;
mkdir(workingDir)
mkdir(workingDir,'movies')
%% create video reader
shuttleVideo = VideoReader('shuttle.avi');
% %%  create image sequence
% ii = 320;
% 
% while hasFrame(shuttleVideo)
%    img = readFrame(shuttleVideo);
%    filename = ['polytest8_',sprintf('%03d',ii),'.jpg'];
%    fullname = fullfile(workingDir,'images',filename);
%    imwrite(img,fullname)    % Write out to a JPEG file (img1.jpg, img2.jpg, etc.)
%    ii = ii+1;
% end
 
%% Find Image File Names
% imageNames = dir(fullfile(workingDir,'images','vort_DiagUV*.png'));
imageNames = dir(fullfile(workingDir,['./images/','vort_DiagUV*.jpg']));
imageNames = {imageNames.name}';
 
%% Create New Video with the Image Sequence
outputVideo = VideoWriter(fullfile(workingDir,'movies',[casename,'_vort.avi']));
outputVideo.FrameRate = shuttleVideo.FrameRate*0.5;
open(outputVideo)
% Loop through the image sequence, load each image, and then write it to the video.
 
for ii = 1:length(imageNames)
	%    img = imread(fullfile(workingDir,'images',imageNames{ii}));
	   img = imread(fullfile([workingDir,'./images/'],imageNames{ii}));
	      writeVideo(outputVideo,img)
end
      %Finalize the video file.
close(outputVideo)
       


