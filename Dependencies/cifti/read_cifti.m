fn='fdt_matrix2_Row_2109.dscalar.nii';
fid = fopen(fn,'r');
header=struct;

str=fgets(fid);
while isempty(strfind(str,'VolumeDimensions'))
    str=fgets(fid);
end

idx= strfind(str,'"');
header.VolumeDimensions=str2num(str(idx(1)+1:idx(2)-1));

fgets(fid);
header.TransformationMatrix(1,:)=str2num(fgets(fid));
header.TransformationMatrix(2,:)=str2num(fgets(fid));
header.TransformationMatrix(3,:)=str2num(fgets(fid));
str=fgets(fid);
header.TransformationMatrix(4,:)=str2num(str(1:strfind(str,'<')-1));

str=fgets(fid);
while isempty(strfind(str,'</MatrixIndicesMap>'))
    str=fgets(fid);
    
    if  ~isempty(strfind(str,'<BrainModel'))
        
        
        idx=strfind(str,'BrainStructure="');
        endidx=strfind(str(idx:end),'"');
        Structure=str(idx+endidx(1):idx+endidx(2)-2);
        
        idx=strfind(str,'Offset="');
        endidx=strfind(str(idx:end),'"');
        str(idx+endidx(1):idx+endidx(2)-2);
        header.BrainModel.(Structure).Offset = str2num(str(idx+endidx(1):idx+endidx(2)-2));
        
        idx=strfind(str,'Count="');
        endidx=strfind(str(idx:end),'"');
        header.BrainModel.(Structure).Count = str2num(str(idx+endidx(1):idx+endidx(2)-2));
        
        idx=strfind(str,'ModelType="');
        endidx=strfind(str(idx:end),'"');
        header.BrainModel.(Structure).ModelType = str(idx+endidx(1):idx+endidx(2)-2);
        
        
        str=fgets(fid);
        t(1,:)=str2num(str(strfind(str,'>')+1:end));
        i=2;
        str = fgetl(fid);
        while isempty(strfind(str,'<'))
            t(i,:)=str2num(str);
            i=i+1;
            str = fgetl(fid);
        end
        
        if ~isequal(size(t,1),header.BrainModel.(Structure).Count)
            disp('Error, Incorrect number of voxels')
        end
        header.BrainModel.(Structure).VoxelCoordinates = t;
        clear t
    end
end
fclose(fid);

vox=[];
Structures=fieldnames(header.BrainModel);
for i=1:numel(Structures)
    vox=[vox;header.BrainModel.(Structures{i}).VoxelCoordinates];
end
header.AllCoordinates=vox;
% fclose(fn);

%% Reconstruct volume 
 [ cifti ] = ciftiopen(fn,'/homes/sparisot/Documents/HCP/workbench/bin_linux64/wb_command');

Conn=zeros(header.VolumeDimensions);
nonZ=find(cifti.cdata~=0);
for i=1:numel(nonZ)
Conn(header.AllCoordinates(nonZ(i),1),header.AllCoordinates(nonZ(i),2),header.AllCoordinates(nonZ(i),3)) = cifti.cdata(nonZ(i));
end
ImViewer(Conn)


figure;vol3d(Conn)
h=vol3d('cdata',Conn);view(3);vol3d(h)

%%

Conn12= zeros(header.VolumeDimensions);
dist= zeros(header.VolumeDimensions);
LVcoor=zeros(numel(LV),3);
for j=1:numel(LV)
Conn=zeros(header.VolumeDimensions);
nonZ=find(fibres.cdata(LV(j),:)~=0);
fib = fibres.cdata(LV(j),:);
LVcoor(j,:)=header.AllCoordinates(find(fib==max(fib),1,'first'),:);
dist(LVcoor(j,1),LVcoor(j,2),LVcoor(j,3))=1;
for i=1:numel(nonZ)
Conn(header.AllCoordinates(nonZ(i),1),header.AllCoordinates(nonZ(i),2),header.AllCoordinates(nonZ(i),3)) = fib(nonZ(i));
end
Conn12=Conn12+Conn;
clear nonZ Conn
% if intersect(j,[20,50,100])%,numel(LV)])
%     figure;h=vol3d('cdata',Conn12);view(3);vol3d(h)
% end
end
D=bwdistsc(dist);
D=rescale(D,0.1,1);

Connperm=permute(Conn12.*D,[2,1,3]);

[x y z]=ind2sub(size(Connperm),find(Connperm>=500));
% figure;plot3(x, y, z, 'k.');
figure;plot(test,Surf)
alpha(0.5)
hold on
plot3(y, x, z, 'k.','Color','y');

Conn2= zeros(header.VolumeDimensions);
LV2coor=zeros(numel(LV2),3);
dist2= zeros(header.VolumeDimensions);
for j=1:numel(LV2)
Conn=zeros(header.VolumeDimensions);
nonZ=find(fibres.cdata(LV2(j),:)~=0);
fib = fibres.cdata(LV2(j),:);
LV2coor(j,:)=header.AllCoordinates(find(fib==max(fib),1,'first'),:);
dist2(LV2coor(j,1),LV2coor(j,2),LV2coor(j,3))=1;
for i=1:numel(nonZ)
Conn(header.AllCoordinates(nonZ(i),1),header.AllCoordinates(nonZ(i),2),header.AllCoordinates(nonZ(i),3)) = fib(nonZ(i));
end
Conn2=Conn2+Conn;
clear nonZ Conn
% if intersect(j,[20,50,100])%,numel(LV2)])
%     figure;h=vol3d('cdata',Conn2);view(3);vol3d(h)
% end
end
D2=bwdistsc(dist2);
D2=rescale(D2,0.1,1);

Connperm=permute(Conn2.*D2,[2,1,3]);

[x y z]=ind2sub(size(Connperm),find(Connperm>=500));
% figure;plot3(x, y, z, 'k.');
figure;plot(test,Surf)
alpha(0.5)
hold on
plot3(y, x, z, 'k.','Color','y');

%%


B=header.TransformationMatrix^(-1);
test=Brain;
% for i=1:32492
% test.vertices(i,:)=test.vertices(i,:)*B(1:3,1:3);
% end
test.vertices=test.vertices*B(1:3,1:3);
test.vertices=bsxfun(@plus,test.vertices,B(1:3,4)');

Connperm=permute(Conn2.*D2,[2,1,3]);

[x y z]=ind2sub(size(Connperm),find(Connperm>=500));
% figure;plot3(x, y, z, 'k.');
figure;plot(test,Surf)
alpha(0.5)
hold on
plot3(y, x, z, 'k.','Color','y');


figure;plot(test)
alpha(0.2)
h=vol3d('cdata',Connperm);view(3);vol3d(h);

%%

cl=39;

Surf.cdata=zeros(32492,1);
Surf.cdata(find(Map1==cl))=1;
cl2=mode(parcelrdm(find(Map1==cl)));
 cl2=24;
Surf.cdata(find(parcelrdm==cl2))=Surf.cdata(find(parcelrdm==cl2))+2;
figure;plot(Brain,Surf)

LV=find(Map1==cl);
LV2=(find(parcelrdm==cl2));
