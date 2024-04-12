% load toolpath
% xn=303;yn=303;r=3;
% this function takes CL points taken from gouge avoidance module
% and finds the optimum sidestep and forward step to keep scallop
% and forward error approximately constant throughout the toolpath
% finally, it gives a zig zag toolpathfinal which can be postprocessed
% xn and yn are number of voxels in x and y directions respectively,
% r is the radius
function [toolpathfinal]=control_scallop(toolpath,xn,yn,r)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- side steps calculation ---------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%divide toolpath final into 2D matrix representing x,y,z grids
x=(reshape(toolpath(:,1),xn,[]));
y=(reshape(toolpath(:,2),xn,[]));
z=(reshape(toolpath(:,3),xn,[]));
xf=zeros(xn,yn); % initialize final x,y,z matrix
yf=zeros(xn,yn);
zf=zeros(xn,yn);
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
last=size(xf,1);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% work on each x layer to get constant scallop CL in x direction
for l=1:size(x,2) %each x layer i.e. side layers
    ip=1;
    a=[x(:,l),y(:,l),z(:,l)]; %take all the side steps (CL)
%     a=a';
    k=2;
    stpt=a(1,:);
    tp(1,:)=stpt;
    ratio=zeros(length(a),1);
    for i=2:length(a) %find scallop between consecutive side voxels
%        if i-ip>1
        a2=a(i,1);a1=stpt(1);
        b2=a(i,3);b1=stpt(3);
%         a21 = a2-a1; b21 = b2-b1;
%         d2 = a21^2+b21^2;
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
        %phai=acosd((r-0.05)/r);
        %phai=10;
%         ratio=b21/a21;
%         if ratio>5.67
%             talbeta=sqrt(2*r*0.05);
%             theta=atan(ratio);
%             d=0.05*((1+1/r)*(1/0.05)*(cos(theta))^2*(1-sqrt(1+4*tan(theta)*talbeta)+2*tan(theta)*talbeta)-1)+0.05;
%         else
            d=r-sqrt(r^2-((a2-a1)^2+(b2-b1)^2)/4);
%         end
%         yindx=a21/0.02;
%         zindx=b21/0.02;
%         dindx=yindx*301+1+zindx;
%         d=dmat(dindx);
        if d>0.05 %if scallop less than the limit check next voxel
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
    %%%% to add different side paths in one matrix matrix dimension has
    %%%% to match hence, diff is used to add extra zeros.
    if i==length(a)
        tp(k,:)=a(i,:);
        k=k+1;
    end
    diff=size(tp,1)-size(xf,2);
    tp=vertcat(tp,zeros(abs(diff),3));
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

% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---------------------- foward step calculation --------------------%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%         b2=a(i,2);b1=stpt(2);
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
% %%%%%%
% % zz=zf;
% % zlayers=zz;
% % for smooth=1:10; %% smoothening in x direction
% %     for j=1:size(zlayers,2)
% %         for i=2:size(zlayers,1)-1
% %             zlayers(i,j)=(zlayers(i-1,j)+zlayers(i+1,j))/2;
% %         end
% %     end
% % end
% % for smooth=1:10; %% smoothening in y direction
% %     for i=1:size(zlayers,1)
% %         for j=2:size(zlayers,2)-1
% %             zlayers(i,j)=(zlayers(i,j-1)+zlayers(i,j+1))/2;
% %         end
% %     end
% % end
% % zf=zlayers;
% %%%%%%
% for i=1:size(xf,2)
%     if mod(i,2)==0
%         xf(:,i)=flipud(xf(:,i));
%         yf(:,i)=flipud(yf(:,i));
%         zf(:,i)=flipud(zf(:,i));
%     end
% end
% toolpathfinal=[reshape(xf,[],1),reshape(yf,[],1),reshape(zf,[],1)];
% %% toolpath generation
% toolpathfinal(~any(toolpathfinal(:,1),2),:)=[];


%------------------------------------------------------------------%

%%%%%%
% post(toolpathfinal,105);
% post(newsurfacepoints,105);
% end
% parfor i=1:10
%     cv(i,:)=toolpath(i,:)
% end
% load toolpath
% xn=303;yn=303;r=3;
% this function takes CL points taken from gouge avoidance module
% and finds the optimum sidestep and forward step to keep scallop
% and forward error approximately constant throughout the toolpath
% finally, it gives a zig zag toolpathfinal which can be postprocessed
% xn and yn are number of voxels in x and y directions respectively,
% r is the radius
% function [toolpathfinal]=control_scallop(toolpath,xn,yn,r)
% % circle circle intersection is calculated by function circcirc1 instead of
% % matlab library function circicrc(p1,r1,p2,r2)
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---------------------- side steps calculation ---------------------%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%divide toolpath final into 2D matrix representing x,y,z grids
% x=reshape(toolpath(:,1),xn,[]);
% y=reshape(toolpath(:,2),xn,[]);
% z=reshape(toolpath(:,3),xn,[]);
% % b=toolpath;
% % [r,c] = size(a);
% % nlay  = 152;
% % matrix   = permute(reshape(b',[c,r/nlay,nlay]),[2,1,3]);
% xf=zeros(xn,yn); % initialize final x,y,z matrix
% yf=zeros(xn,yn);
% zf=zeros(xn,yn);
% nn=1;
% %%work on each x layer to get constant scallop CL in x direction
% for l=1:size(x,1) %each x layer i.e. side layers
%     a=[x(l,:)',y(l,:)',z(l,:)']; %take all the side steps (CL)
%     k=2;
%     stpt=a(1,:);
%     tp(1,:)=stpt;
%     for i=2:length(a) %find scallop between consecutive side voxels
%         [xout,yout] = circcirc1(stpt(2),stpt(3),r,a(i,2),a(i,3),r);
%         if stpt(3)==a(i,3)
%             [v,n]=min(yout);
%             pt=[a(i,1),xout(n),v];
%         else
%             [v,n]=min(xout);
%             pt=[a(i,1),v,yout(n)];
%         end
% %         [v,n]=min(yout); %two circle intersection points
%         %%%%%% find tangent
%         a2=a(i,2);a1=stpt(2);
%         b2=a(i,3);b1=stpt(3);
%         a21 = a2-a1; b21 = b2-b1;
%         d2 = a21^2+b21^2;
%         r1=r;
%         r2=r1;
%         r21 = (r2-r1)/d2;
%         s21 = sqrt(d2-(r2-r1)^2)/d2; % <-- If d2<(r2-r1)^2, no solution is possible
%         u1 = [-a21*r21-b21*s21,-b21*r21+a21*s21]; % Left unit vector
%         u2 = [-a21*r21+b21*s21,-b21*r21-a21*s21]; % Right unit vector
%         L1 = [a1,b1]+r1*u1; L2 = [a2,b2]+r2*u1; % Left line tangency points
%         R1 = [a1,b1]+r1*u2; R2 = [a2,b2]+r2*u2; % Right line tangency points
%         %%%%%%min distance from tangent to the intersection point (scallop)
%         %         d = point_to_line([a(i,1),xout(n),v], [a(i,1),R1(1),R1(2)], [a(i,1),R2(1),R2(2)]);
%         if L1(1,1)<R1(1,1)
%             d = point_to_line(pt, [a(i,1),L1(1),L1(2)], [a(i,1),L2(1),L2(2)]);
%         else
%             d = point_to_line(pt, [a(i,1),R1(1),R1(2)], [a(i,1),R2(1),R2(2)]);
%         end
%         if d<0.05 || isnan(d)==1%if scallop less than the limit check next voxel
%             nn=nn+1;
%             continue
%         else %else save the voxel as sidestep
%             tp(k,:)=a(i-1,:);
%             ss(k)=(nn);
%             k=k+1;
%             nn=1;
%             stpt=a(i-1,:);%and calculate side step from this new voxel
%         end
%     end
%     %%%% to add different side paths in one matrix matrix dimension has
%     %%%% to match hence, diff is used to add extra zeros.
%     diff=length(tp)-size(xf,2);
%     tp=vertcat(tp,zeros(abs(diff),3));
%     xf(l,:)=tp(:,1)';
%     yf(l,:)=tp(:,2)';
%     zf(l,:)=tp(:,3)';
%     tp=[];
%     a=[];
% end
% %%remove all zeros columns form the final x,y,z matrix
% xff=yf;
% xf( :, ~any(xff,1) ) = [];  %columns
% yf( :, ~any(xff,1) ) = [];  %columns
% zf( :, ~any(xff,1) ) = [];  %columns
% %%%%%% smoothening
% % zz=zf;
% % zlayers=zz;
% % for smooth=1:10; %% smoothening in x direction
% %     for j=1:size(zlayers,2)
% %     for i=2:size(zlayers,1)-1
% %         zlayers(i,j)=(zlayers(i-1,j)+zlayers(i+1,j))/2;
% %     end
% %     end
% % end
% % for smooth=1:10; %% smoothening in y direction
% %     for i=1:size(zlayers,1)
% %     for j=2:size(zlayers,2)-1
% %         zlayers(i,j)=(zlayers(i,j-1)+zlayers(i,j+1))/2;
% %     end
% %     end
% % end
% % zf=zlayers;
% %%%%%% flip consectutive paths upside down to perform zig-zag machining
% % for i=1:size(xf,2)
% %     if mod(i,2)==0
% %     xf(:,i)=flipud(xf(:,i));
% %     yf(:,i)=flipud(yf(:,i));
% %     zf(:,i)=flipud(zf(:,i));
% %     end
% % end
% toolpathfinal=[reshape(xf,[],1),reshape(yf,[],1),reshape(zf,[],1)];
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %---------------------- foward step calculation --------------------%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%         [xout,yout] = circcirc1(stpt(1),stpt(2),r,a(i,1),a(i,2),r);
%         if stpt(1)==a(i,1) %if consecutive voxels lie on same x then minimum of x is found
%             [v,n]=min(xout);
%             pt=[v,a(i,3),yout(n)]; %yout gives the z coordinate
%         else
%             [v,n]=min(yout); % if x is different for consecutive voxels then find min of z
%             pt=[xout(n),a(i,3),v];
%         end
% %         [v,n]=min(yout);
%         %%%%%% find tangent
%         a2=a(i,1);a1=stpt(1);
%         b2=a(i,2);b1=stpt(2);
%         a21 = a2-a1; b21 = b2-b1;
%         d2 = a21^2+b21^2;
%         r1=r;
%         r2=r1;
%         r21 = (r2-r1)/d2;
%         s21 = sqrt(d2-(r2-r1)^2)/d2; % <-- If d2<(r2-r1)^2, no solution is possible
%         u1 = [-a21*r21-b21*s21,-b21*r21+a21*s21]; % Left unit vector
%         u2 = [-a21*r21+b21*s21,-b21*r21-a21*s21]; % Right unit vector
%         L1 = [a1,b1]+r1*u1; L2 = [a2,b2]+r2*u1; % Left line tangency points
%         R1 = [a1,b1]+r1*u2; R2 = [a2,b2]+r2*u2; % Right line tangency points
%         %%%%%%min distance from tangent to the intersection point
%         %         d = point_to_line([xout(n),a(i,2),v], [R1(1),a(i,2),R1(2)], [R2(1),a(i,2),R2(2)]);
%         if L1(1,2)<R1(1,2)
%             d = point_to_line(pt, [L1(1),a(i,3),L1(2)], [L2(1),a(i,3),L2(2)]);
%         else
%             d = point_to_line(pt, [R1(1),a(i,3),R1(2)], [R2(1),a(i,3),R2(2)]);
%         end
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
%     diff=length(tp)-size(xf,1);
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
% % zz=zf;
% % zlayers=zz;
% % for smooth=1:10; %% smoothening in x direction
% %     for j=1:size(zlayers,2)
% %         for i=2:size(zlayers,1)-1
% %             zlayers(i,j)=(zlayers(i-1,j)+zlayers(i+1,j))/2;
% %         end
% %     end
% % end
% % for smooth=1:10; %% smoothening in y direction
% %     for i=1:size(zlayers,1)
% %         for j=2:size(zlayers,2)-1
% %             zlayers(i,j)=(zlayers(i,j-1)+zlayers(i,j+1))/2;
% %         end
% %     end
% % end
% % zf=zlayers;
% %%%%%%
% for i=1:size(xf,2)
%     if mod(i,2)==0
%         xf(:,i)=flipud(xf(:,i));
%         yf(:,i)=flipud(yf(:,i));
%         zf(:,i)=flipud(zf(:,i));
%     end
% end
% toolpathfinal=[reshape(xf,[],1),reshape(yf,[],1),reshape(zf,[],1)];
% %% toolpath generation
% toolpathfinal(~any(toolpathfinal(:,1),2),:)=[];
% %%%%%%
% % post(toolpathfinal,105);
% % post(newsurfacepoints,105);
% % end
% % parfor i=1:10
% %     cv(i,:)=toolpath(i,:)
% % end
