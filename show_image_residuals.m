
function show_image_residuals(f, measured, estimated);
    global vidObj1;
    figure(f);
    cla(gca);
    len = length(measured);
    tmp1 = reshape(measured, 2, len/2)';
    tmp2 = reshape(estimated, 2, len/2)';
    xlim([4000, 6000]);
    ylim([4000, 6000]);
    plot(tmp1(:,1), tmp1(:,2), 'gx');
    hold on;
    plot(tmp2(:,1), tmp2(:,2), 'bx');
    legend('measured projections', 'estimated projections');   
  %  writeVideo(vidObj1, getframe(f));
end