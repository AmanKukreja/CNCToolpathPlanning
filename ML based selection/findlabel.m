function [modellabelADA1,modellabelISO1,modellabelHYB1]=findlabel(px,py,pz,n,vs,r)
% reshape x, y, and z 2d matrices to column matrix and pass surface points
% to the toolpath generation codes

% create a steep feature
x =round(reshape(px,[],1),2);
y =round(reshape(py,[],1),2);
z =round(reshape(pz,[],1),2)+10;

surfacepts =[x y z];

[toolpathfinal]=adaptive_planar(surfacepts,n,n,r);
toolpathfinal(toolpathfinal(:,3)==0,:)=[];
[modellab]=labels(toolpathfinal,r,vs);
modellabelADA1=modellab;

% subplot(2,2,2)
% plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))

[toolpathfinal]=control_scallop(surfacepts,n,n,r);
toolpathfinal(toolpathfinal(:,3)==0,:)=[];
[modellab]=labels(toolpathfinal,r,vs);
modellabelISO1=modellab;

% subplot(2,2,3)
% plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))

try
    [MCP,RecReg,RecReg1]=Dynamic_clustering_for_MCP(surfacepts,n,vs);
    if isempty(MCP)==1
        modellabelHYB1=zeros(1,6);
    else
        [toolpathfinal,surfaceptsregion]=control_scallop_reg_st(surfacepts,n,n,r,RecReg,RecReg1);
        toolpathfinal(toolpathfinal(:,3)==0,:)=[];

%         subplot(2,2,4)
%         plot3(toolpathfinal(:,1),toolpathfinal(:,2),toolpathfinal(:,3))
%         hold on
        
        maxz=max(toolpathfinal(:,3));
        regionpath=[toolpathfinal(size(toolpathfinal,1),1),toolpathfinal(size(toolpathfinal,1),2),maxz];
        % plot regions toolpath
        for i=1:length(surfaceptsregion)
            toolpathregion=surfaceptsregion{i};
            toolpathregion(toolpathregion(:,3)==0,:)=[];            
            regionpath=vertcat(regionpath,[regionpath(size(regionpath,1),1),regionpath(size(regionpath,1),2),maxz],[toolpathregion(1,1),toolpathregion(1,2),maxz],toolpathregion);

%             plot3(toolpathregion(:,1),toolpathregion(:,2),toolpathregion(:,3))
            %             hold on
        end
        toolpathfinal=vertcat(toolpathfinal,regionpath);

        [modellab]=labels(toolpathfinal,r,vs);
        modellabelHYB1=modellab;
%         hold off
    end
catch
    modellabelHYB1=zeros(1,6);
end
