function plotting(w,h,wd,hd,w1,h1,degree)
plot(w,h,'LineWidth',2);
hold on
plot(wd,hd,'LineWidth',2);
plot(w1,h1,'LineWidth',2);
legend('Estimated','Desired','Equiripple')
title(['Degree= ' num2str(degree)])
drawnow;
hold off
end

