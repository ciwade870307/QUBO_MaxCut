function [] = plot_acceptProb(arg,mode_gif)

figure('position',[400,300,1000,300]); hold on; box on; grid on;
set(gca,'FontSize',12); set(gca,'looseInset',[0 0 0 0]);

x_axis = 0:arg.maxStep-1;
h2 = plot(x_axis, arg.p_list,'.black', 'MarkerSize',10);
h2 = plot(x_axis, arg.p_list,'.black', 'MarkerSize',10); % Duplicate

xlabel("Step")
ylabel("Acceptance Probability")
title("Acceptance Probability");

F_name = "DA@had3_step_"+arg.maxStep;
filename = "../Fig/qap/"+F_name+".gif";

if mode_gif
    t_interval = 0.05;
    for idx_step = 1:arg.maxStep
        set(h2, 'xdata',idx_step-1, 'ydata', arg.p_list(idx_step),'MarkerSize',20, 'color','red');		% 設定新的 y 座標
        title("Acceptance Probability, step="+idx_step);
        pause(t_interval)
        drawnow

    %     Save as fig
        frame = getframe(gcf);
        img =  frame2im(frame);
        [img,cmap] = rgb2ind(img,256);
        if idx_step == 1
            imwrite(img,cmap,filename,'gif','LoopCount',Inf,'DelayTime',t_interval);
        else
            imwrite(img,cmap,filename,'gif','WriteMode','append','DelayTime',t_interval);
        end
    end
end

end