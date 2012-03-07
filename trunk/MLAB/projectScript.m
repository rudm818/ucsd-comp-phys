%% simple test for collapsing sphere (project code outline)
%
% AUTHOR: Michael Folkerts
% EMAIL: mfolkerts@physics.ucsd.edu
% WEBSITE: http://bit.ly/folkerts
% DATE: Mar. 2012
%
% This scrip:
% 1) Generates a uniform spherical distribution
% 2) Saves the data as an Aarseth input file
% 3) Call's the Aarseth code and specifies the input file
% 4) Then calls my previous script to read the output and create an animated plot.
% 

%_PARAMETERS_
workingDirectory = './'; % where my scripts are
% the following paths should be relative to the 'workingDirectory'
inputFilePath =  '../data/sphere.txt';
outputFilePath = 'aarseth.data'; % this was defined in Aarseth code
pathToAarseth =  '../bin/aarseth';

% for plot/animation
axisMax = 1;

% Aarseth Parameters (remember G = 1)
numberOfPts =    1000;
massOfPts =      0.5; % in your units
eta =            0.02; % ~ precision
timeStep =       0.005; % in your units
finalTime =      0.5;
epsilonSquared = 0.25;


%% _START_SCRIPT_
% make sure we are in the working folder
cd(workingDirectory)

% generate iniform sphere data
% NOTE: one could read in galaxy data here instead
maxRadius = .5;
rShift = [0,.5,.5]; % shift for all points
vCM = [0,-1.75,-1.75]; % drift for all points

R = maxRadius*rand(numberOfPts,1);
theta = pi*rand(numberOfPts,1);
phi = 2*pi*rand(numberOfPts,1);

% .* is element wise multiplication (not dot product)
x = R.*cos(phi).*sin(theta);
y = R.*sin(phi).*sin(theta);
z = R.*cos(theta);

r = horzcat(x,y,z); % group and transpose to make a list of row vectors
              % ex: r(1,:) = position (row) vector of the first particle 
              %     r(1,2) = the y value of position for first particle 
% set zero velocity for center of mass:
v = zeros(numberOfPts,3); % list of row vectors


% Transform velocity and position
% (there might be a more elegant way to do this)
% (one could also rotate a coordinate system as well) 
for i = 1:numberOfPts
    r(i,:) = r(i,:) + rShift; % shift all points by vector rShift
    v(i,:) = v(i,:) + vCM; % add velocity of center of mass to each velocity vector
end

% TODO: plot to verify data here

% save Aarseth input file:
fileID = fopen(inputFilePath,'w');

% Header: numParticles, eta, timeStep, finalTime and epsilonSquared
fprintf(fileID,'%i %f %f %f %f\n',numberOfPts, eta, timeStep, finalTime, epsilonSquared); 

for i = 1:numberOfPts
    fprintf(fileID,'%f %f %f %f %f %f %f\n',massOfPts,r(i,1),r(i,2),r(i,3),v(i,1),v(i,2),v(i,3));
end

fclose(fileID);

% clean up memory?
% clear r v theta phi x y z

%% call Aarseth code
system([pathToAarseth ' < ' inputFilePath],'-echo');

%%
% call my previous scrip (as a function)
animate3D(outputFilePath,numberOfPts,axisMax)
