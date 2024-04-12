function [fitresult, gof, scallop, ii] = createFit13(xx, yy, ii, scallop,x,y)
%Curve fitting over center points of intersection % create HISTOGRAM for the whole surface using distance
%approach
%% Fit: 'Fittedcurve'.
[xData, yData] = prepareCurveData( x, y ); %fitting on the center points
% [xData, yData] = prepareCurveData( xx, yy ); %fitting on the center points
% Set up fittype and options. curve fitting according to the number of intersetions
% if length(x)<10
%     fits=['poly' num2str(length(x)-1)];
% else
%     fits='poly9';
% end
% ft = fittype( fits );
ft = fittype( 'smoothingspline' );
% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft,'Normalize', 'on' );
%% minimum distance to curve calculation
kk=1;
curvexy=[xx fitresult(xx)];
mapxy=[xx yy];
[~,distance] = distance2curve(curvexy,mapxy,'linear');
% plot(xData,distance);
%% convert distance below fitted curve to negative
% for i=1:length(distance)
%     if output.residuals(i)<0
%         distance(i)=-1*distance(i);
%     end
% end
%% find minimum lowest points between scallops
for sc=50:length(distance)-50
    if distance(sc)<distance(sc-1)&&distance(sc)<distance(sc-2)&&distance(sc)<distance(sc+1)&&distance(sc)<distance(sc+2)
        ptofint(kk)=sc;
        kk=kk+1;
    end
end
%% Scallop calculation
for i=1:kk-2 %find intervals for scallop calculation based on side step
    scallop(ii)= max(distance(ptofint(i):ptofint(i+1)))-min(distance(ptofint(i):ptofint(i+1)));
    %sidesteps(ii)= xx(ptofint(i+1))-xx(ptofint(i));
    ii=ii+1;
end


