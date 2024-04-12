% surfacepts(:,1:2)=surfacepts(:,1:2)*0.02; %either convert all into double or all into integer
% clear all
% clc
% load('0.02_Face_Alper_full.mat')
function [MCP,RecReg,RecReg1]=Dynamic_clustering_for_MCP(surfacepts,n,vs)
% for i=1:length(surfacepts)
%     if surfacepts(i,2)~=surfacepts(i+1,2)
%         break
%     end
% end
% xd=i;
% yd=length(surfacepts)/xd;
xd=n; yd=n; % taken as the model we are considering is a square model.. to generate more general we will have to change this.
r=3;
%% reduce dimenions of the voxel grid 50 times for 1 mm step size
z=reshape(surfacepts(:,3),xd,[]);
y=reshape(surfacepts(:,2),xd,[]);
x=reshape(surfacepts(:,1),xd,[]);

k=1;
for i=1:10:size(z,1)
    z1(k,:)=z(i,:);
    k=k+1;
end
k=1;
for i=1:10:size(z,1)
    y1(k,:)=y(i,:);
    k=k+1;
end

k=1;
for i=1:10:size(z,1)
    x1(k,:)=x(i,:);
    k=k+1;
end

k=1;
for i=1:10:size(z,2)
    z2(:,k)=z1(:,i);
    y2(:,k)=y1(:,i);
    x2(:,k)=x1(:,i);
    k=k+1;
end

%% Remove voxels that are outside boundary from getting selected as a region
% for irem=1:size(z2,1)
%     k=1;
%     for jrem=1:size(z2,2)
%         if (sum(z2(irem,1:jrem))==0)
%             savejrem(k)=jrem;
%             value=z2(irem,jrem+1);
%             k=k+1;
%         end
%     end
% %     [~,index]=min(nonzeros(z2(irem,:)));
%     z2(irem,savejrem)=value;
%     savejrem=[];
%     k=1;
%     for jrem=1:size(z2,2)
%         if (sum(z2(irem,jrem:size(z2,2)))==0)
%             value=z2(irem,jrem-1);
%             break
%         end
%     end
%     z2(irem,jrem:size(z2,2))=value;
% end
% z2(z2==0)=300;
% % comment 78, uncomment 79 to select highest grad region
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find gradient of the surface in x and y directions
[gx,~]=gradient(z2);
gx=abs(atand(gx)); % here gx is along y axis as during reshaping the matrix gets transformed.
if max(gx,[],'all')<45 %if gradient more that 45 degrees then no segmentation
    MCP=[];
    RecReg=[];
    RecReg1=[];
    return
end
% we are taking gradient along y because it is side step direction. If we
% want to generate toolpath along x then we should make it gy.
%% find the cluster of various gradients using dynamic k-means techinque
[lb,center] = adaptcluster_kmeans(gx);
% lb(lb==length(center))=1;
%Note: we are ignoring the largest gradient region because in non-rectangular components boundary points give high gradient

%% Find the clusters within the maximum gradient cluster, this will
%aid us in finding out the MCP
im=y2;
% im(lb~=length(center)-1)=0;
im(lb~=length(center))=0;
[lb,center] = adaptcluster_kmeans(im);
% surf(lb)
%% Plot the surface and final cluster centers, it signifies the main gradient regions
X=reshape(x2,[],1);
Y=reshape(y2,[],1);
Z=reshape(z2,[],1);

% plot3(X,Y,Z,'.');
% surf(z1)
maxz=max(max(z));
% hold on
k=1;

for i=2:length(center) % first cluster is of gradient 0 that is why ignored
    regy=y2(lb==i);
    centery=mean(regy);
    regx=x2(lb==i);
    centerx=mean(regx);
    xl=[min(regx) max(regx)];
    yl=[min(regy) max(regy)];

    if i==2 || (yl(2)-RecReg{k-1}(2,2))>5 || ((centerx-centerxo)>5)
        RecReg1{k}=[round((xl-x(1,1)+vs)/vs);round((yl-y(1,1)+vs)/vs)];
        RecReg{k}=[xl;yl];
        %     plot3(centerx,centery,maxz,'ro')
        centerxo=centerx;
        k=k+1;
    else
        xl=[min(xl(1),RecReg{k-1}(1,1)) max(xl(2),RecReg{k-1}(1,2))];
        yl=[min(yl(1),RecReg{k-1}(2,1)) max(yl(2),RecReg{k-1}(2,2))];
        RecReg1{k-1}=[round((xl-x(1,1)+vs)/vs);round((yl-y(1,1)+vs)/vs)];
        RecReg{k-1}=[xl;yl];
        centerxo=centerx;
    end
    x1=[];
    y1=[];
end

% plot3(xl,yl,[100,100],'bo')

MCP=[round((centerx-x(1,1)+vs)/vs),round((centery-y(1,1)+vs)/vs)];
% % % % For Adaptive planar band % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RecReg=[round((xl-x(1,1)+0.02)/0.02);round((yl-y(1,1)+0.02)/0.02)];
% MCP=[x(round((centerx-x(1,1)+0.02)/0.02),1),y(1,round((centery-y(1,1)+0.02)/0.02))];

end
%% remove unwanted boundary 0s from face alper model
% mz=min(min(z));
% for i=1:size(z,2)
%     for j=1:500
%         if z(j,i)==0
%             z(j,i)=mz;
%         end
%     end
% end
% for i=1:size(z,2)
%     for j=6000:6525
%         if z(j,i)==0
%             z(j,i)=mz;
%         end
%     end
% end
%
% k=1;
% for i=1:50:size(z,1)
%     z1(k,:)=z(i,:);
%     k=k+1;
% end
%
% k=1;
% for i=1:50:size(z,2)
%     z2(:,k)=z1(:,i);
%     k=k+1;
% end