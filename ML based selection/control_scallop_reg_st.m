% This is for face model right now, it removes the high gradient region
% points and stitch the remaining points around the region to generate the
% iso-scallop toolpath throughout. then adaptive planar toolpath can be
% generated in the regions separately. This code simply lifts the tool
% up when at high gradient region
function [toolpathfinal,surfaceptsregion]=control_scallop_reg_st(toolpath,xn,yn,r,RecReg,RecReg1)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- side steps calculation ---------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%divide toolpath final into 2D matrix representing x,y,z grids
x=(reshape(toolpath(:,1),xn,[]));
y=(reshape(toolpath(:,2),xn,[]));
z=(reshape(toolpath(:,3),xn,[]));
nn=1;
d=0;
error=0;

% xs=0:0.02:6;
% ys=0:0.02:6;
% comb=zeros(90601,2);
% ks=1;
% for i=1:length(xs)
%     for j=1:length(ys)
%         comb(ks,:)=[xs(i) ys(j)];
%         ks=ks+1;
%     end
% end
% dmat=r-sqrt(r^2-(comb(:,1).^2+comb(:,2).^2)/4);

%% if you want to keep constant forward step use this code, this speeds up
%%%% the computation also as forward step is less significant.
k=1;
last=size(x,1);
for i=1:10:last
    xfn(k,:)=x(i,:);
    yfn(k,:)=y(i,:);
    zfn(k,:)=z(i,:);
    k=k+1;
end
xfn(k,:)=x(last,:);
yfn(k,:)=y(last,:);
zfn(k,:)=z(last,:);

x=xfn';
y=yfn';
z=zfn';

% xf=zeros(size(x)); % initialize final x,y,z matrix
% yf=zeros(size(x));
% zf=zeros(size(x));
xf=zeros(xn,yn); % initialize final x,y,z matrix
yf=zeros(xn,yn);
zf=zeros(xn,yn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% work on each x layer to get constant scallop CL in x direction
for l=1:size(x,2) %each x layer i.e. side layers
    ip=1;
    a=[x(:,l),y(:,l),z(:,l)]; %take all the side steps (CL)
    %     a=a';
    k=2;
    stpt=a(1,:);
    tp(1,:)=stpt;
    %     ratio=zeros(length(a),1);
    for i=2:length(a) %find scallop between consecutive side voxels
        %        if i-ip>1
        a2=a(i,1);a1=stpt(1);
        b2=a(i,3);b1=stpt(3);
        a21 = a2-a1; b21 = b2-b1;
        d2 = a21^2+b21^2;
        % Comment it if you want to create toolpath without considering intermediate voxels
        %         if d>0.04 %correct scallop for points having scallop near to
        %         desired value by considering intermediate voxels
        %             xin=a(ip:i,2);
        %             yin=a(ip:i,3);
        %             m=-b21/a21;
        %             c=m*a1-b1;
        %             error=max(abs(m*xin+yin+c)/norm([m,1]));
        %             ratio=(yin-stpt(3))./(xin-stpt(2));
        %             [maxrat,indrat]=max(ratio);
        %             ratio(i)=b21/a21;
        %             [maxrat,indrat]=max(ratio(ip+1:i-1));
        %             theta=atan(maxrat);
        %             thetas=atan(ratio(i));
        %             maxd=sin(theta-thetas);
        %             li=(a(indrat+ip,2)-stpt(2));
        %             error=maxd*li/cos(theta);
        %             if abs(error)>0.05
        %                 d=0.06;
        %                 error=0;
        %             else
        %                 d=r-sqrt(r^2-d2/4)+error;
        %             end
        %        end
        %         else
        %             d=r-sqrt(r^2-d2/4);
        %         end
        %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        d=r-sqrt(r^2-d2/4);
        %         yindx=a21/0.02;
        %         zindx=b21/0.02;
        %         dindx=yindx*301+1+zindx;
        %         d=dmat(dindx);
   
        
        if d>=0.05 %if scallop less than the limit check next voxel
            if round(d,3)>0.05
                tp(k,:)=a(i-1,:);
            else
                tp(k,:)=a(i,:);
            end
            d=0;
            %             ss(k)=(nn);
            k=k+1;
            %             nn=1;
            stpt=tp(k-1,:);%and calculate side step from this new voxel
            ip=i;
        else %else save the voxel as sidestep
            %             nn=nn+1;
            continue
        end
    end
    %%%% to add different side paths in one matrix, matrix dimension has
    %%%% to match. hence, diff is used to add extra zeros.
    if i==length(a)
        tp(k,:)=a(i,:);
        k=k+1;
    end
    diff=size(tp,1)-size(xf,2);
    tp=vertcat(tp,zeros(abs(diff),3));
%     if size(xf,2)<size(tp,1)
%         tp=tp(1:size(xf,2),:);
%     else
%         tp=vertcat(tp,zeros(abs(diff),3));
%     end
    xf(l,:)=tp(:,1)';
    yf(l,:)=tp(:,2)';
    zf(l,:)=tp(:,3)';
    tp=[];
    a=[];
end
%%remove all zeros columns form the final x,y,z matrix
xff=yf;
xf( :, ~any(xff,1) ) = [];  %columns
yf( :, ~any(xff,1) ) = [];  %columns
zf( :, ~any(xff,1) ) = [];  %columns

xff=yf;
xf( ~any(xff,2),: ) = [];  %columns
yf( ~any(xff,2),: ) = [];  %columns
zf( ~any(xff,2),: ) = [];  %columns

%% Region segmentation and toolpath generation
surfaceptsregion=[];
% xf(44:92,19:43)=0;
for i=1:length(RecReg)
    while true
        [~,col]=find((yf>RecReg{i}(2,1))&(xf>RecReg{i}(1,1))); % select the start point of the high gradient region
        col=min(col);
        if col==[]
            RecReg{i}(2,1)=RecReg{i}(2,1)-vs; % it is an iterative process, as the point may not be present in the existing array.
        else
            startregy=max(2,abs(col)-2); % extra paths taken in y(side step direction) to include curved paths of isoscsllop
            break % here 5 side steps are taken
        end
    end
    col=[];
    while true
        [~,col]=find(yf>RecReg{i}(2,2)&(xf>RecReg{i}(1,1))); % end point column
        col=min(col);
        if col==[]
            RecReg{i}(2,2)=RecReg{i}(2,2)+vs;
        else
            endregy=min(col+1,size(yf,2)); % steps ahead of the regions
            break
        end
    end
    startregx=max(2,abs(find(xf==RecReg{i}(1,1),1,'first')-2)); % taking 5 steps extra in start and end points in fwd direction
    endregx=min(size(x,2),(find(xf==RecReg{i}(1,2),1,'first')+2));
    
    yff=yf;% to pick scallop points which went inside the region and are to be used for next layers
    zff=zf;
    xff=xf;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yf(startregx:endregx,startregy:endregy)=0; %remove all the points in the region
    zf(startregx:endregx,startregy:endregy)=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Generate Adaptive Planar for Region %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    startregyr=find(y(:,1)==min(yf(:,startregy-1)),1,'first'); %Find the start side path by finding out the closest point on iso-scallop toolpath
    endregyr=find(y(:,1)==max(yf(:,endregy+1)),1,'first'); % from the base matrix x,y,z on which toolpath will be generated
    
    xr=x(startregyr:endregyr,startregx:endregx)'; % extract points of region on which toolpath is to be generated
    yr=y(startregyr:endregyr,startregx:endregx)';
    zr=z(startregyr:endregyr,startregx:endregx)';
    
    %     xr=x(RecReg1{i}(2,1)-500:RecReg1{i}(2,2)+500,startregx:endregx)';
    %     yr=y(RecReg1{i}(2,1)-500:RecReg1{i}(2,2)+500,startregx:endregx)';
    %     zr=z(RecReg1{i}(2,1)-500:RecReg1{i}(2,2)+500,startregx:endregx)';
    regionpoints=[reshape(xr,[],1), reshape(yr,[],1), reshape(zr,[],1)];
    if isempty(xr)
        break
    end
    [toolpathReg{i}]=adaptive_planar_w(regionpoints,size(xr,1),size(xr,2),r);
    surfaceptsregion=toolpathReg{i};
    divv=endregx-startregx+1; % number of side paths.
    yrr=reshape(surfaceptsregion(:,2),divv,[]); % convert toolpath back to matrix form
    xrr=reshape(surfaceptsregion(:,1),divv,[]);
    zrr=reshape(surfaceptsregion(:,3),divv,[]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %     for ii=startregx:endregx
    %         yf(ii,startregy:endregy)=yf(startregx-1,startregy:endregy); %keeping their x points same add lift to the region points
    %         zf(ii,startregy:endregy)=30;
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%% Find Overlapping Paths, stitch them to the Iso-scallop
    %%%%%%%%%%%%%% toolpath and remove it from adaptive planar patch path
    for ii=startregy:endregy
        j=2;
        while j<size(yrr,2)
            if (yrr(1,j)>=yf(startregx-1,ii))&&(yrr(1,j-1)<=yf(startregx-1,ii))
                if abs(yrr(1,j)-yf(startregx-1,ii))>abs(yrr(1,j-1)-yf(startregx-1,ii))
                    indexr=j-1;
                    break
                else 
                    indexr=j;
                    break
                end
            else
                j=j+1;
            end
        end
%         if indexr>size(yrr,2)
            yf(startregx:endregx,ii)=yrr(:,min(indexr,size(yrr,2))); %keeping their x points same add lift to the region points
            zf(startregx:endregx,ii)=zrr(:,min(indexr,size(yrr,2)));
            yrr(:,min(indexr,size(yrr,2)))=[];
            xrr(:,min(indexr,size(yrr,2)))=[];
            zrr(:,min(indexr,size(yrr,2)))=[];
%         end
    end
    
    %% slide the iso-scallop toolpath points ahead of region, close to the region boundary
    start=abs(find(y(:,1)==yf(startregx,endregy-1),1,'first')); % taking 5 steps extra in start and end points in fwd direction
    points=[reshape(x(start:yn,startregx:endregx)',[],1), reshape(y(start:yn,startregx:endregx)',[],1), reshape(z(start:yn,startregx:endregx)',[],1)];
    ydd=yn-start+1;xdd=endregx-startregx+1;
    [poix,poiy,poiz]=control_scallop_1(points,xdd,ydd,r);
    sz=size(yf,2)-endregy+2-size(poix,2);sx=endregx-startregx+1;
    poix=horzcat(poix,zeros(sx,sz));poiy=horzcat(poiy,zeros(sx,sz));poiz=horzcat(poiz,zeros(sx,sz));
    xf(startregx:endregx,endregy-1:size(yf,2))=poix;
    yf(startregx:endregx,endregy-1:size(yf,2))=poiy;
    zf(startregx:endregx,endregy-1:size(yf,2))=poiz;
    
%     for ii=startregx:endregx %we want to delete the points and drag the limit to the region boundary
%         for j=endregy-5:size(yf,2)
%             if yff(ii,j)<yf(startregx-1,endregy+1)
%                 continue
%             else
%                 if (yff(ii,j)-yf(startregx-1,endregy+1))<abs((yff(ii,j-1)-yf(startregx-1,endregy+1)))
%                     jk=j;
%                 else
%                     jk=j-1;
%                 end
%                 tempy=yff(ii,jk:size(yff,2));
%                 yf(ii,endregy+1:size(yf,2))=0;
%                 yf(ii,endregy+1:size(yff,2)-(jk-(endregy+1)))=tempy;
%                 tempx=xff(ii,jk:size(xff,2));
%                 xf(ii,endregy+1:size(xf,2))=0;
%                 xf(ii,endregy+1:size(xff,2)-(jk-(endregy+1)))=tempx;
%                 tempz=zff(ii,jk:size(zff,2));
%                 zf(ii,endregy+1:size(zf,2))=0;
%                 zf(ii,endregy+1:size(zff,2)-(jk-(endregy+1)))=tempz;
%                 break
%             end
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Regpoints{i}=[{xrr}, {yrr}, {zrr}];
    xrr=[];yrr=[];zrr=[];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% yf(44:92,19:43)=0;
% zf(44:92,19:43)=0;
%
% for i=44:92
%     yf(i,19:43)=yf(43,19:43);
%     zf(i,19:43)=30;
% end
%
% for i=44:92
%     for j=44:127
%         if yf(i,j)<97
%             continue
%         else
%             tempy=yf(i,j:127);
%             yf(i,43:127)=0;
%             yf(i,44:127-(j-44))=tempy;
%             tempx=xf(i,j:127);
%             xf(i,43:127)=0;
%             xf(i,44:127-(j-44))=tempx;
%             tempz=zf(i,j:127);
%             zf(i,43:127)=0;
%             zf(i,44:127-(j-44))=tempz;
%             break
%         end
%     end
% end
%%
%%%%%% smoothening
% zz=zf;
% zlayers=zz;
% for smooth=1:10; %% smoothening in x direction
%     for j=1:size(zlayers,2)
%     for i=2:size(zlayers,1)-1
%         zlayers(i,j)=(zlayers(i-1,j)+zlayers(i+1,j))/2;
%     end
%     end
% end
% for smooth=1:10; %% smoothening in y direction
%     for i=1:size(zlayers,1)
%     for j=2:size(zlayers,2)-1
%         zlayers(i,j)=(zlayers(i,j-1)+zlayers(i,j+1))/2;
%     end
%     end
% end
% zf=zlayers;

% Abruptboundaries can be removed here
% abt=10; %for face model
% xf=xf(abt:size(xf,1)-abt,:);
% yf=yf(abt:size(yf,1)-abt,:);
% zf=zf(abt:size(zf,1)-abt,:);

%%%%%% flip consectutive paths upside down to perform zig-zag machining
for i=1:size(xf,2)
    if mod(i,2)==0
        xf(:,i)=flipud(xf(:,i));
        yf(:,i)=flipud(yf(:,i));
        zf(:,i)=flipud(zf(:,i));
    end
end
clear surfaceptsregion
%%%%%% Flip consectutive paths upside down to perform zig-zag machining for
%%%%%% region points
if exist('Regpoints','var')
    for reg=1:length(Regpoints)
        for i=1:size(Regpoints{reg}{1},2)
            if mod(i,2)==0
                xrr(:,i)=flipud(Regpoints{reg}{1}(:,i));
                yrr(:,i)=flipud(Regpoints{reg}{2}(:,i));
                zrr(:,i)=flipud(Regpoints{reg}{3}(:,i));
            else
                xrr(:,i)=Regpoints{reg}{1}(:,i);
                yrr(:,i)=Regpoints{reg}{2}(:,i);
                zrr(:,i)=Regpoints{reg}{3}(:,i);
            end
        end
        surfaceptsregion{reg}=[reshape(xrr,[],1) reshape(yrr,[],1) reshape(zrr,[],1)];
        xrr=[];yrr=[];zrr=[];
    end
else
    surfaceptsregion=0;
end
toolpathfinal=[reshape(xf,[],1),reshape(yf,[],1),reshape(zf,[],1)];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------- foward step calculation according to adaptive planar--------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clearvars -except toolpathfinal xn yn r
% x=reshape(toolpathfinal(:,1),xn,[]); %toolpath again reshaped to x,y,z
% y=reshape(toolpathfinal(:,2),xn,[]);
% z=reshape(toolpathfinal(:,3),xn,[]);
% yj=size(x,2);%final maximum sidestep given by the no. of columns of x
% xf=zeros(xn,yj);
% yf=zeros(xn,yj);
% zf=zeros(xn,yj);
% nn=1;
% %%work on each y layer to get constant forward error toolpath
% for l=1:size(x,2) % move in forward direction
%     a=[x(:,l),y(:,l),z(:,l)];
%     k=2;
%     stpt=a(1,:);
%     tp(1,:)=stpt;
%     for i=2:length(a)
%         a2=a(i,1);a1=stpt(1);
%         b2=a(i,3);b1=stpt(3);
%         a21 = a2-a1; b21 = b2-b1;
%         d2 = a21^2+b21^2;
%         d=r-sqrt(r^2-d2/4);
%         if d<0.05 || isnan(d)%forward error limit
%             nn=nn+1;
%             continue
%         else
%             tp(k,:)=a(i-1,:);
%             ss(k)=(nn);
%             k=k+1;
%             nn=1;
%             stpt=a(i-1,:);
%         end
%     end
%     if i==length(a)
%         tp(k,:)=a(i,:);
%         k=k+1;
%     end
%     diff=size(tp,1)-size(xf,1);
%     tp=vertcat(tp,zeros(abs(diff),3));
%     xf(:,l)=tp(:,1);
%     yf(:,l)=tp(:,2);
%     zf(:,l)=tp(:,3);
%     tp=[];
%     a=[];
% end
% xff=xf;
% xf( ~any(xff,2),: ) = [];  %rows
% yf( ~any(xff,2),: ) = [];  %rows
% zf( ~any(xff,2),: ) = [];  %rows
% %%%%%%
% for i=1:size(xf,2)
%     if mod(i,2)==0
%         xf(:,i)=flipud(xf(:,i));
%         yf(:,i)=flipud(yf(:,i));
%         zf(:,i)=flipud(zf(:,i));
%     end
% end
% toolpathfinal=[reshape(xf,[],1),reshape(yf,[],1),reshape(zf,[],1)];
%% toolpath generation
toolpathfinal(~any(toolpathfinal(:,1),2),:)=[];