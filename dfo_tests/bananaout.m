function stop = bananaout(x,~,state)
persistent ph
stop = false;
switch state
    case 'init'
       newplot
       xx = -2:.2:2;
       yy = -1:.2:3;
       [xx,yy] = meshgrid(xx,yy);
       zz = 100*(yy-xx.^2).^2+(1-xx).^2;
       % Set up the appropriate colormap
       % In this case, the colormap has been chosen to give the surf plot
       % a nice healthy banana color.
       hsv2 = hsv;
       hsv3 = [hsv2(11:64,:);hsv2(1:10,:)];
       % draw the surf plot
       surface(xx,yy,zz,'EdgeColor',[0.8,0.8,0.8]);
       xlabel 'x(1)'
       ylabel 'x(2)'
       view(10,55);
       colormap(hsv3);
       hold on;
       [~,contHndl] = contour3(xx,yy,zz,[100,500],'k');
       contHndl.Color = [0.8,0.8,0.8];
       
       plot3(-1.9,2,267.62,'ko','MarkerSize',15,'LineWidth',2);
       text(-2,2.4,267.62,'   Start','Color',[0 0 0]);
       plot3(1,1,0,'ko','MarkerSize',15,'LineWidth',2);
       text(0.8,1.4,0,'   Solution','Color',[0 0 0]);
       drawnow

    case 'iter'
        x1 = x(1);
        y1 = x(2);
        z1 = 100*(y1-x1.^2).^2+(1-x1).^2;
        ph = plot3(x1,y1,z1,'r.','MarkerSize',25);
        h = gca;
        h.SortMethod = 'childorder';
        drawnow;
        
    case 'done'
        legend(ph,'Iterative steps','Location','east')
        hold off
end