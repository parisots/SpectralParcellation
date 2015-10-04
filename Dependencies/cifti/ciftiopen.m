function [ cifti ] = ciftiopen(filename,caret7command)
%Open a CIFTI file by converting to GIFTI external binary first and then
%using the GIFTI toolbox

if strfind(filename,'wbsparse')
    idx=strfind(filename,'.');
    unix([caret7command  ' -convert-matrix4-to-matrix2 ' filename ' ' filename(1:idx) 'dconn.nii']);
    filename = [filename(1:idx) 'dconn.nii'];
end

grot=fileparts(filename);
if (size(grot,1)==0)
grot='.';
end
tmpname = tempname(grot);

tic
unix([caret7command ' -cifti-convert -to-gifti-ext ' filename ' ' tmpname '.gii']);
toc

tic
cifti = gifti([tmpname '.gii']);
toc

unix(['rm ' tmpname '.gii ' tmpname '.gii.data']);

end

