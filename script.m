clear all
close all

K = [1,1]+rand(1,2);
K = K ./ norm(K) * 9;


gamma = 5;
gamma_sq = gamma^2;

circ = 0:0.05:2*pi+0.06;
x_circ = [cos(circ); sin(circ)];

f = figure('Position', [250, 250, 1000, 1200]);
p1 = polarplot(circ, gamma*(circ.*0+1), 'r:', 'linewidth', 3);
ax = gca();
ax.RLim = [0 10];
ax.RTick = [0:10];
hold on
p2 = polarplot(ax, circ, gamma*(circ.*0+1), 'b', 'linewidth', 4);
legend({'Circle R=5', 'K*x sweep'}, 'fontsize', 15)
grid on
t1 = annotation('textbox', [0.70, 0.23, 0.1, 0.1], 'String', 't',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);
t2 = annotation('textbox', [0.70, 0.19, 0.1, 0.1], 'String', 't',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);
t3 = annotation('textbox', [0.70, 0.33, 0.1, 0.1], 'String', 't',...
'FaceAlpha',.9, 'backgroundcolor', 'w', 'interpreter', 'latex', 'fontsize', 15);
t4 = annotation('textbox', [0.70, 0.29, 0.1, 0.1], 'String', 't',...
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

% t1 = text(pi*0.12, 16, 't1');
% t2 = text(pi*0.2, 17, 't1');
% t3 = text(pi*0.25, 20, 't1');
% t4 = text(pi*0.28, 23, 't1');
%%
k = 100;
frames = [];
while k > 3

    K = round(K.*0.99, 2);
    lm = [1, K;
        K', eye(2)*gamma_sq];

    qm = eye(2)*gamma_sq - K'*K;

    lmeigs = round(eigs(lm), 2);
    qmeigs = round(eigs(qm), 2);
    Keig = round(sqrt(eigs(K'*K)), 2);  
    
    Kx = K*x_circ;

        
    p2.RData = abs(Kx);
    t1.String = strcat("$eigs(Schur) : ", strjoin(string(lmeigs), ' || ' ), '$');
    t2.String = strcat("$eigs(\gamma^2 I - K^TK): ", strjoin(string(qmeigs), ' || ' ), '$');
    t3.String = strcat("$\sqrt{eigs(K^TK)} : ", strjoin(string(Keig), ' || ' ), '$');
    t4.String = strcat("$K : ", strjoin(string(K), ' || ' ), '$');
    
    
    k = max(Kx);
    frames = [frames, getframe(f)];
end

writerObj = VideoWriter('animation');
writerObj.FrameRate = 5;

open(writerObj);
for k = 2:length(frames)
   frame = frames(k);
   writeVideo(writerObj, frame);
end
close(writerObj);