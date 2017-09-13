close all;

nTfs = length(TF.cell_list);
numPoints = 50;

for jTf=1:nTfs
    tf = TF.cell_list{jTf};
    fun = tf.fun;
    %fun = TF.f_addNoiseBig(tf.fun);
    if ( tf.dim == 2 || tf.dim == Inf )
        fig = figure('name',tf.name);
        [lb,ub,xbst] = TF.f_getVectors(tf, 2);
        x = lb(1) : (ub(1)-lb(1))/numPoints : ub(1);
        y = lb(2) : (ub(2)-lb(2))/numPoints : ub(2);
        [x,y] = meshgrid(x,y);
        z = zeros(size(x,1),size(x,2));
        for j=1:size(x,1)
            for k=1:size(x,2)
                z(j,k)=fun([x(j,k);y(j,k)]);
            end
        end
        surfc(x,y,z,'EdgeColor',[0.8,0.8,0.8]);
        alpha(0.5);
        xlabel('x');
        ylabel('y');
        zlabel('fval');
        hold on;

        plot3(xbst(1),xbst(2),tf.fbst,'ko','MarkerSize',15,'LineWidth',2);

        saveas(fig, [pwd '/images/' tf.name '.png']); 
    end
end