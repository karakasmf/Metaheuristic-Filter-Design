function plotting(w,h,wd,hd,degree)
plot(w,h,'LineWidth',2);
hold on
plot(wd,hd,'LineWidth',2);
legend('Estimated','Desired')
title(['Degree= ' num2str(degree)])
drawnow;
hold off
end

