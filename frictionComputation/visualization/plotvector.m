% function plotCOR(alpha,beta,x0,y0,z0)
function plotvector(a,b,c,x0,y0,z0,color,linewidth)

if nargin<8
    linewidth = 1.5;
end
%     a = cos(alpha) * cos(beta);
%     b = cos(alpha) * sin(beta);
%     c = sin(alpha);

    direction = [a;b;c];
    location = [x0;y0;z0];

    P0CORGlobal = location - 0.5 * direction;
    P1CORGlobal = location + 0.5 * direction;
    vectarrow(P0CORGlobal,P1CORGlobal,color,linewidth);
    hold on
    plot3(x0,y0,z0,'r.','MarkerSize',5);


%     plot3(P(1,:),P(2,:),P(3,:),'b.');  
end
