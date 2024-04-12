function [modellabels]=labels(toolpathfinal,r,vs)
NCS=[toolpathfinal(:,2), toolpathfinal(:,1), toolpathfinal(:,3)];
ii=1;
scallop=[];
%% Toolpath length
TPL=0;
for i=1:length(NCS)-1
    TPL=TPL+norm(NCS(i,:)-NCS(i+1,:));
end
TPL=TPL/1000;

% xlabel('X-axis (mm)');
% ylabel('Y-axis (mm)');
% zlabel('Z-axis (mm)');
% set(findall(gcf,'-property','FontSize'),'FontSize',18)
% center=[(max(NCS(:,1))+min(NCS(:,1)))/2 (max(NCS(:,2))+min(NCS(:,2)))/2];
% NCS=[NCS(:,1:2)-repmat(center,length(NCS),1) NCS(:,3)];
% NCS=NCS(10:length(NCS)-10,:);
% NCS=NCS*R;
% plot3(NCS(:,1),NCS(:,2),NCS(:,3),'b')
xmax=max(toolpathfinal(:,1));
xmin=min(toolpathfinal(:,1));
xlen=xmax-xmin+vs;
ylen=max(toolpathfinal(:,2))-min(toolpathfinal(:,2))+vs;
value =xmin:(xlen/40):xmax ; %number of planes

%% Plane intersections
for jj=2:length(value)-5
    n=[1 0 0]; %normal in the diection of x
    V0=[value(jj) 0 0]; %plane at point = value
    kk=1;
    for i=1:length(NCS)-1
        P0=NCS(i,:);
        P1=NCS(i+1,:);
        [I,check]=plane_line_intersect(n,V0,P0,P1); %check specifies if the plane intersects at a unique point inside P0P1 or not
        if check==1%(I(1,1)<P0(1,1) && I(1,1)>P1(1,1)) || (I(1,1)>P0(1,1) && I(1,1)<P1(1,1))
            inters(kk,:)=I;
            kk=kk+1;
        end
    end
    kk=1;
    if exist('inters')==0
        continue
    end
    %     plot3(inters(:,1),inters(:,2),inters(:,3))
    %     inters(length(inters)-10:length(inters),:)=[];
    %     inters(inters(:,3)==0,:)=[];
    %     inters=sortrows(inters,2);
    %     plot3(inters(:,1),inters(:,2),inters(:,3))
    sizein=size(inters);
    if sizein(1,1)<3
        continue
    end
    rt=r; %RADIUS OF THE TOOL
    for i=1:sizein(1,1)-1
        [xout,yout] = circcirc(inters(i,2),inters(i,3),rt,inters(i+1,2),inters(i+1,3),rt);
        if yout(1,2)<yout(1,1) %choose lower intersection point
            sp(kk,:)=[xout(1,2),yout(1,2)]; %intersection scallop points
        else
            sp(kk,:)=[xout(1,1),yout(1,1)];
        end
        kk=kk+1;
    end
    %     plot(sp(:,1),sp(:,2),'.')
    %%generation of scallop curves
    xx=[];
    yy=[];
    for i=2:sizein(1,1)-1
        z=inters(i,3);
        y=inters(i,2);
        z1=sp(i-1,2);
        y1=sp(i-1,1);
        z2=sp(i,2);
        y2=sp(i,1);
        t1=90+atand((y1-y)/(z-z1)); %initial angle for one scallop
        if y<y2
            t2=180-atand((z-z2)/(y2-y)); %final angle
        else
            t2=90-atand((y-y2)/(z-z2));
        end
        xc=(y-rt*cosd(t1:.9:t2))'; % scallop points found by forming arc at side steps points
        yc=(z-rt*sind(t1:.9:t2))';
        %plot(xc,yc,'.')
        xx=vertcat(xx,xc); %merge previous scallop points with the next side step scallop points
        yy=vertcat(yy,yc);
        xc=[];yc=[];
    end
    A=[xx yy];
    A = A( ~any( isnan( A ) | isinf( A ), 2 ),: ); % remove inf and NAN points
    xx=A(:,1);
    yy=A(:,2);
    x=inters(:,2);
    y=inters(:,3);
    %           [fitresult, gof, scallop, ii,sidesteps] = createFit1(xx, yy, ii, scallop,sidesteps); %curve fitting for all scallops at once
    [fitresult, gof, scallop, ii] = createFit13(xx, yy, ii, scallop,x,y);
    %         [fitresult, gof, scallop] = createFit23(x,y,xx, yy);
    %         [fitresult, gof, scallop] = createFit2(xx, yy);
    % scallop=[];
    %clear fitresult gof output A inters sp line
    xx=[];
    yy=[];
    A=[];
    inters=[];
end
%% Results
% r=3;
% vs=0.02;
% xlen=max(NCS(:,1))-min(NCS(:,1))+vs;
% ylen=max(NCS(:,2))-min(NCS(:,2))+vs;
h=0.05;
% scallop(scallop<0.03)=[];
scallop(scallop>0.08)=[];
greater=sum(scallop>h)/length(scallop);
avg=h-abs(mean(scallop)-h);
deviation=std(scallop);
% maximum side step is taken as for a flat surface
ms=sqrt(4*(2*r*h-h^2));
TPLscore=(xlen*ylen/(ms*1000))/TPL;
Abruptness=0;

for i=1:length(NCS)-2
    if round(NCS(i+2,1),1)~=round(NCS(i,1),2)
        if ((NCS(i+2,1)==NCS(i,1))&&(NCS(i+2,2)==NCS(i,2)))||((NCS(i+1,1)==NCS(i,1))&&(NCS(i+1,2)==NCS(i,2)))
            continue
        else
            t2=atan((NCS(i+2,2)-NCS(i,2))/(NCS(i+2,1)-NCS(i,1)));
            t1=atan((NCS(i+1,2)-NCS(i,2))/(NCS(i+1,1)-NCS(i,1)));
            t3=(NCS(i+2,:)-NCS(i,:));
            d3=sqrt(t3(1)^2+t3(2)^2);
            Abruptness=Abruptness+abs((NCS(i+1,1)-NCS(i,1))*(sin(t1-t2)/cos(t1)))/d3;
        end
    end
end

% for i=1:length(NCS)-2
% %     if round(NCS(i+2,1),1)~=round(NCS(i,1),2) %x(3) not equal to x(1)
%         if ((NCS(i+2,1)==NCS(i+1,1))&&(NCS(i+2,2)~=NCS(i+1,2)))||((NCS(i+1,1)==NCS(i,1))&&(NCS(i+1,2)~=NCS(i,2)))
%             Abruptness=Abruptness+0.5;
%             continue %x(2)==x(1)&&y(2)==y(1)
%         else
%             t1=(NCS(i+2,:)-NCS(i,:));
%             d1=sqrt(t1(1)^2+t1(2)^2);
%             t2=(NCS(i+2,:)-NCS(i+1,:));
%             d2=sqrt(t2(1)^2+t2(2)^2);
%             t3=(NCS(i+1,:)-NCS(i,:));
%             d3=sqrt(t3(1)^2+t3(2)^2);
% %             d1/(d3+d2)
% %             plot(NCS(i:i+2,1),NCS(i:i+2,2),'ro')
%             if (d2+d3)==0
%                 Abruptness=Abruptness+1;
%             else
%                 Abruptness=Abruptness+d1/(d3+d2);
%             end
%         end
% %     end
% end

Abruptness=1-Abruptness/(length(NCS)-2); % sir method A/arc length
% Abruptness=Abruptness/(length(NCS)-2); % our new proposed method base/arc_length
% metric to specify a strategy as better than the other one
% TotalScore=0.3*TPLscore+0.3*avg/h-0.1*(greater)-deviation-0.3*Abruptness;
TotalScore=0.3333*TPLscore+0.3333*avg/h+0.3333*Abruptness;
TotalScore1=0.3*TPLscore+0.3*avg/h+0.4*Abruptness;
% modellabels=[TotalScore,TPLscore,avg/h,greater,deviation,Abruptness];
modellabels=[TotalScore,TotalScore1,TPLscore,avg/h ,Abruptness];
% TotalScore=0.3*a(1)+0.3*a(2)-0.1*(a(3))-a(4)-0.3*a(5);