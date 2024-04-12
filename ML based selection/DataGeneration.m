% %% convert voxel based surface model into 2d array and save as a csv file
% for i=1:length(surfacepts)
%     if surfacepts(i,2)~=surfacepts(i+1,2)
%         break
%     end
% end
% xd=i;
% sp=1/0.02;
% yd=length(surfacepts)/xd;
% x=reshape(surfacepts(:,3),xd,yd);
% x1=x(1:sp:xd,:);
% x2=x1(:,1:sp:yd);
% filename=input("Enter File Name",'s');
% writematrix(x2,[filename,'.csv'])
% 
% %% convert data generated in 1000X1000 to 100X100 to pass onto training
% indexxes1=1:10:1000;
% for i=1:1000
%     temp=check(indexxes1,:,i);
%     modelfortraining(:,:,i)=temp(:,indexxes1);
% end
% datafortraining=reshape(modelfortraining,[],1000);
% 
% writematrix(finalmodelt,'dataset_a.csv')
% writematrix(finallabels,'labels_a.csv')
% % datafortraining_shaped=reshape(datafortraining,100,100,3000);
% 
% %% %%%%%%%%%%%%% Data Visualization code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% indexxes=10:10:1000000;
% indexxes1=10:10:1000;
% hybridlabels=find(modelmetrics(:,3)~=0);
% % x=reshape(x,[],1);
% % y=reshape(y,[],1);
% for kk=1:length(hybridlabels)
%     i=hybridlabels(kk);
% %     surface=reshape(model(:,kk),1000,[])';
%     surface=model(indexxes,i);
%     surface=reshape(surface,100,[]);
%     surface=surface(:,indexxes1)';
%     subplot(2,2,1)
%     surf(surface)
%     if labelofbest(i)==1
%         best='Preffered: Adaptive';
%     elseif labelofbest(i)==2
%         best='Preffered: Iso-scallop';
%     else
%         best='Preffered: Hybrid';
%     end
% 
%     title(best);
% 
%     surfacepts=[x,y,model(:,i)];
%     [toolpathfinal]=adaptive_planar(surfacepts,n,n,r);
%     toolpathfinal(toolpathfinal(:,3)==0,:)=[];
%     subplot(2,2,2)
%     plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))
%     title(directionADA(i))
% 
%     [toolpathfinal]=control_scallop(surfacepts,n,n,r);
%     toolpathfinal(toolpathfinal(:,3)==0,:)=[];
%     subplot(2,2,3)
%     plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))
%     title(directionISO(i))
%     if modelmetrics(i,3)==0
%         subplot(2,2,4)
%         plot3(0,0,0)
%         continue
%     end
%     [MCP,RecReg,RecReg1]=Dynamic_clustering_for_MCP(surfacepts,n,vs);
%     try
%         if isempty(MCP)==1
%             modellabelsHYB(k,:)=zeros(1,6);
%         else
%             [toolpathfinal,surfaceptsregion]=control_scallop_reg_st(surfacepts,xd,yd,r,RecReg,RecReg1);
%             toolpathfinal(toolpathfinal(:,3)==0,:)=[];
%             subplot(2,2,4)
%             plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))
%             hold on
%             maxz=max(toolpathfinal(:,3));
%             regionpath=[toolpathfinal(size(toolpathfinal,1),1),toolpathfinal(size(toolpathfinal,1),2),maxz];
%             % plot regions toolpath
%             for i=1:length(surfaceptsregion)
%                 toolpathregion=surfaceptsregion{i};
%                 toolpathregion(toolpathregion(:,3)==0,:)=[];
%                 plot3(toolpathregion(:,1),toolpathregion(:,2),toolpathregion(:,3))
%             end
%         end
%     catch
%         modellabelHYB1=zeros(1,6);
%     end
%     null=[];
% end
% % Plot Only hybrid models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hybridlabels=find(modelmetrics(:,3)~=0);
% % x=reshape(x,[],1);
% % y=reshape(y,[],1);
% for kk=1:length(hybridlabels)
%     i=hybridlabels(kk);
%     surface=reshape(model(:,kk),100,[])';
%     surf(surface)
%     if labelofbest(i)==1
%         best=['Preffered: Adaptive ' num2str(directionADA(i)) ];
%     elseif labelofbest(i)==2
%         best=[ 'Preffered: Iso-scallop ' num2str(directionISO(i))];
%     else
%         best=['Preffered: Hybrid ', num2str(directionHYB(i))];
% 
%     end
% 
%     title(best);
% 
% end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Data generation code %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc; clear all; clf;
% close all;
indexxes1=1:10:1000;
vs=0.1;
r=3;
check=zeros(36,1);
%p = xlsread('bzdata_mane2.xlsx');  % check this input surface in region boundaryfn also
p=zeros(36,3);
load('sixpointszeros.mat');
px=zeros(1000,1000);
py=px;
pz=px;
model=zeros(10000,1000);
modellabelsADA=zeros(1000,6);
directionADA=zeros(1000,1);
modellabelsISO=zeros(1000,6);
directionISO=zeros(1000,1);
modellabelsHYB=zeros(1000,6);
directionHYB=zeros(1000,1);
f=waitbar(0,'Creating Dataset');
for k=1:10000
    waitbar(k/10000,f)
%     try
        numr=randi([1,3]); % number of peeks
        for i=1:numr
            pl=randi([1,36]); % placement of control points of peaks
            p(pl,3)=randi([0,100]); % height of control points of peak
        end
        check(:,k)=p(:,3); % save control points of all the models
        siz=6; %number of control points in one direction
        cpx=reshape(p(:,1),siz,[])';
        cpy=reshape(p(:,2),siz,[])';
        cpz=reshape(p(:,3),siz,[])';

        % code for generation of basis functions of bezier surface of degree n
        % n=siz-1;
        % syms u v f(u)
        % for i=1:n+1
        %     k=i-1;
        %     bu(i)=(factorial(n)/(factorial(k)*factorial(n-k)))*u^k*(1-u)^(n-k);
        %     bv(i)=(factorial(n)/(factorial(k)*factorial(n-k)))*v^k*(1-v)^(n-k);
        % end
        % f(u)=bu;
        % f(v)=bv;

        if mod(k,1000)==0
            pause(120)
        end

        ux = 0.001:0.001:1;
        vx = 0.001:0.001:1;
        [m,n] = size(ux);
        for i = 1:n
            for j = 1:n
                u = ux(i);
                v = vx(j);
                bu = [-(u - 1)^5 5*u*(u - 1)^4 -10*u^2*(u - 1)^3 10*u^3*(u - 1)^2 -5*u^4*(u - 1) u^5];
                bv = [-(v - 1)^5 5*v*(v - 1)^4 -10*v^2*(v - 1)^3 10*v^3*(v - 1)^2 -5*v^4*(v - 1) v^5];
                px(i,j)= bu*cpx*bv'; % actual surface points in x y and z
                py(i,j)= bu*cpy*bv';
                pz(i,j) = bu*cpz*bv';
            end
        end
        % surf(px,py,pz);
        %x=px;y=py;z=pz;

             
        temp=pz(indexxes1,indexxes1);
        zm =round(reshape(temp,[],1),2)+10;
        model(:,k)=zm; % save models
        
        % toolpath towards positive y direction [1]
        [modellabelADA1,modellabelISO1,modellabelHYB1]=findlabel(px,py,pz,n,vs,r);
        z=rot90(pz); % toolpath towards positive x direction [2]
        [modellabelADA2,modellabelISO2,modellabelHYB2]=findlabel(px,py,z,n,vs,r);
        z=rot90(z); % toolpath towards negative y direction [3]
        [modellabelADA3,modellabelISO3,modellabelHYB3]=findlabel(px,py,z,n,vs,r);
        z=rot90(z); % toolpath towards negative x direction [4]
        [modellabelADA4,modellabelISO4,modellabelHYB4]=findlabel(px,py,z,n,vs,r);
        [modellabelsADA(k,:),directionADA(k)]=max([modellabelADA1(1),modellabelADA2(1),modellabelADA3(1),modellabelADA4(1)]);
        [modellabelsISO(k,:),directionISO(k)]=max([modellabelISO1(1),modellabelISO2(1),modellabelISO3(1),modellabelISO4(1)]);
        [modellabelsHYB(k,:),directionHYB(k)]=max([modellabelHYB1(1),modellabelHYB2(1),modellabelHYB3(1),modellabelHYB4(1)]);
         
        %     filename=['model' num2str(k)];
        %     save(filename,'points')
        p(:,3)=0;
%     catch
%         p(:,3)=0;
%         continue
%     end
    % plot3(surfacepts(:,1),surfacepts(:,2),surfacepts(:,3))
    % plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))
    % for Iso-scallop
    % [toolpathfinal]=control_scallop_binary(surfacepts,xd,yd,r); % this code also converts the surface points into zig-zag toolpath
    % toolpathfinal(toolpathfinal(:,3)==0,:)=[];
    % plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))
end
% plot3(points(:,1),points(:,2),points(:,3),'.');
modelmetrics=[modellabelsADA(:,1),modellabelsISO(:,1),modellabelsHYB(:,1)];
[~,labelofbest]=max(modelmetrics,[],2);
sum(labelofbest==1)
sum(labelofbest==2)
sum(labelofbest==3)
save('set10.mat','check', 'directionADA', 'directionISO', 'directionHYB', 'labelofbest', 'model', 'modellabelsADA', 'modellabelsISO', 'modellabelsHYB', 'modelmetrics')