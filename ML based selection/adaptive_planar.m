%%%%%% Adaptive planar without considering intermediate voxels for scallop
%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% It also calculate forward steps

% this function takes CL points (here toolpath) taken from gouge avoidance module
% and finds the optimum sidestep and forward step to generate adptive planar toolpath zig-zag toolpath
% xn and yn are number of voxels in x and y directions respectively,
% r is the radius

%%%% It works well for simple rectangular parts %%%%%%%%%%%%%%%%%%%%%%%%%%%

function [toolpathfinal]=adaptive_planar(toolpath,xn,yn,r)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- side steps calculation ---------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%divide toolpath final 2D matrix representing x,y,z grids
x=reshape(toolpath(:,1),xn,[]);
y=reshape(toolpath(:,2),xn,[]);
z=reshape(toolpath(:,3),xn,[]);
tic

%%%%%%%%%%%%%%%%
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

x=xfn;
y=yfn;
z=zfn;

% xf=zeros(size(x)); % initialize final x,y,z matrix
% yf=zeros(size(x));
% zf=zeros(size(x));
% xf=zeros(xn,yn); % initialize final x,y,z matrix
% yf=zeros(xn,yn);
% zf=zeros(xn,yn);

%%%%%%%%%%%%%%%%
xf(:,1)=x(:,1);
yf(:,1)=y(:,1);
zf(:,1)=z(:,1);
stpt=[x(:,1),y(:,1),z(:,1)]; %start point
k=2;l=2;lp=2;counter=1;
%%work on each x layer to get constant scallop CL in x direction
while l<=size(x,2) %each x layer i.e. side layers
    a=[x(:,l),y(:,l),z(:,l)]; %take all the side steps (CL)
    for i=1:length(a) %find scallop between consecutive side voxels
        
        a2=a(i,1);a1=stpt(i,1);
        b2=a(i,3);%-r;
        b1=stpt(i,3);%-r;
        a21 = a2-a1; b21 = b2-b1;
        if (b1==0)||(b2==0) % if one of b1 or b2 is 0 that means no CC point hence no scallop -> d2=0 -> d=0
            d2=0; %as d2 becomes more than r, scallop becomes imaginary hence, it is corrected to 0
        else
            %c21=a(i,2)-stpt(i,2);
            d2 = (a21^2+b21^2);d2=d2/4;
        end
        
        if b21>=2
            d=r-sqrt(r^2-d2)-(r-a21); % to incorporate slanted features we approximate the value of scallop
        else
            d=r-sqrt(r^2-d2); % d = scallop height
        end
        
        if d<0.05 || isnan(d)%if scallop less than the limit check next voxel
            if i==length(a)
                l=l+1; % increment if all the points on side step have scallop less the desired
            else
                continue
            end
        else  %else save the voxel as sidestep
            %%%%%%%% to ignore small portions of abrupt surfaces counter is
            %%%%%%%% used, here 100 equivalent to a total of 2 mm of abruptness
            %%%%%%%% in the surface is ignored to get the next side step. also,
            %%%%%%%% minimum spacing is kept at 25 steps(0.5mm) to reduce
            %%%%%%%% computations.
%             if counter<100
%                 counter=counter+1;
%                 continue
%             end
%             counter=1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if l==lp % if imidiately next path gives greater scallop then change the index to save next path, basically there wont be any segmentation
                l=l+1;
            end
            xf(:,k)=x(:,l-1);
            yf(:,k)=y(:,l-1);
            zf(:,k)=z(:,l-1);
            stpt=[xf(:,k),yf(:,k),zf(:,k)];
            k=k+1;
            lp=l;
            break
        end
    end
    if l==size(x,1) %for the top most point only
        xf(:,k)=x(:,l);
        yf(:,k)=y(:,l);
        zf(:,k)=z(:,l);
        l=l+1;
    end
    %%%% to add different side paths in one matrix matrix dimension has
    %%%% to match hence, diff is used to add extra zeros.
    %     diff=length(tp)-size(xf,2);
    %     tp=vertcat(tp,zeros(abs(diff),3));
    %     xf(l,:)=tp(:,1)';
    %     yf(l,:)=tp(:,2)';
    %     zf(l,:)=tp(:,3)';
    %     tp=[];
    %     a=[];
end
toc
toolpathfinal=[reshape(xf,[],1),reshape(yf,[],1),reshape(zf,[],1)];
%%%%%%%%%% yf= z coordinate of the toolpath, zf= y coordinate
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- foward step calculation --------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic
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
%         if d<0.05 || isnan(d)==1%forward error limit
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
% toc
% %%%%%%
for i=1:size(xf,2)
    if mod(i,2)==0
        xf(:,i)=flipud(xf(:,i));
        yf(:,i)=flipud(yf(:,i));
        zf(:,i)=flipud(zf(:,i));
    end
end
toolpathfinal=[reshape(xf,[],1),reshape(yf,[],1),reshape(zf,[],1)];
%% toolpath generation
toolpathfinal(~any(toolpathfinal(:,1),2),:)=[];
