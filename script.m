clear all
close all

% generate a random row vector ([1,1] ensures that none of the K gain is 0)
K = [1,1]+rand(1,2);
% normalize to a unit vector, then make the norm 9
K = K ./ norm(K) * 9;

% gamma = max u=Kx assuming |x|_2=1
gamma = 5;
gamma_sq = gamma^2;

% vector of angles 0 to 2pi
circ = 0:0.05:2*pi+0.06;

f = figure('Position', [250, 250, 1000, 1200]);

% Plot a circle with radius gamma
% gamma*(circ.*0 + 1) creates a vector of gammas with circ's size
p1 = polarplot(circ, gamma*(circ.*0+1), 'r:', 'linewidth', 3);

ax = gca();

ax.RLim = [0 10];
ax.RTick = [0:10];
hold on
grid on

% Set up the plot for Kx sweep curve
p2 = polarplot(ax, circ, gamma*(circ.*0+1), 'b', 'linewidth', 4);

legend({'Circle R=5', 'K*x sweep'}, 'fontsize', 15)


% Initialize textboxes
t1 = annotation('textbox', [0.70, 0.23, 0.1, 0.1], 'String', 'init',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);
t2 = annotation('textbox', [0.70, 0.19, 0.1, 0.1], 'String', 'init',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);
t3 = annotation('textbox', [0.70, 0.33, 0.1, 0.1], 'String', 'init',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);
t4 = annotation('textbox', [0.70, 0.29, 0.1, 0.1], 'String', 'init',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);
t5 = annotation('textbox', [0.70, 0.37, 0.1, 0.1], 'String', '$\gamma = 5$',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);

str = "Schur LMI : $ \pmatrix{1 & K \cr K^T & \gamma^2 I} > \mathbf{0}$";
t6 = annotation('textbox', [0.60, 0.80, 0.1, 0.1], 'String', str,...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);

str = "QMI constraint : $ \gamma^2 I - K^T K > \mathbf{0} $";
t7 = annotation('textbox', [0.60, 0.72, 0.1, 0.1], 'String', str,...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);

t7 = annotation('textbox', [0.80, 0.5, 0.1, 0.1], 'String', '$|x|_2=1$',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);

%% Animation
k = 100;
frames = [];
% while max value of Kx is >3
while k > 3

    % decrease K norm
    K = round(K.*0.99, 2);
    
    % Schur LMI
    lm = [1, K;
        K', eye(2)*gamma_sq];

    % QMI
    qm = eye(2)*gamma_sq - K'*K;

    % Calculate eigenvalues
    lmeigs = round(eigs(lm), 2);
    qmeigs = round(eigs(qm), 2);
    Keig = round(sqrt(eigs(K'*K)), 2);  
    
    % Calculate Kx sweep curve
    Kx = K*x_circ;

    % Update the Kx sweep curve plot
    p2.RData = abs(Kx);
    
    % Update the textboxes
    t1.String = strcat("$eigs(Schur) : ", strjoin(string(lmeigs), ' || ' ), '$');
    t2.String = strcat("$eigs(\gamma^2 I - K^TK): ", strjoin(string(qmeigs), ' || ' ), '$');
    t3.String = strcat("$\sqrt{eigs(K^TK)} : ", strjoin(string(Keig), ' || ' ), '$');
    t4.String = strcat("$K : ", strjoin(string(K), ' || ' ), '$');
    
    
    k = max(Kx);
    frames = [frames, getframe(f)];
end

%% Save the animation
writerObj = VideoWriter('animation');
writerObj.FrameRate = 5;

open(writerObj);
for k = 2:length(frames)
   frame = frames(k);
   writeVideo(writerObj, frame);
end
close(writerObj);
