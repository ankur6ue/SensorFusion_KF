function draw_rotating_cube(pos, rot, f)
% for making a video
global vidObj1;
% clear current plot
figure(f);
cla(gca);

% Define a six row by four column matrix to define the six cube faces
fm = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];

% Define an eight row by three column matrix to define the vertices at which
% the faces meet
scale = 0.3;
vm = scale*[0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
vm = vm + pos';

% Plot the cube ----- gives each face a different color and creates the
% cube at a convenient viewing angle
clear cdata;
cdata =     [0 0 0; % black
            1 0 0; % red
            1 0 1; % magenta
            0 0 1; % blue
            0 1 0; % green
            1 1 0; % yellow
            1 1 1; % white
            0 1 1; % cyan
            ];
        
vm = vm*rot;
patch('Vertices',vm,'Faces',fm,'FaceVertexCData',cdata,'FaceColor','interp');

axis equal;
% Limits have to be set based on the size of the desired movement box
xlim(2*[-0.5 0.5]); ylim(2*[-0.5 0.5]); zlim(2*[-0.5 0.5]);
xlabel('x'); ylabel('y'); zlabel('z');
view(3)
writeVideo(vidObj1, getframe(f));
end