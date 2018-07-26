% Author: Veronica Corona, vc324@cam.ac.uk

% Compute convex segmentation and segmentation error. 
% lambda is the regularisation parameter
% k is the weight for prior knowledge

function [segmentation,u,error_segmentation,confusion_matrix] = convex_segmentation(lambda,groundtruth,k,post)

%load post3d;
[~,u,ud]=segmentation(lambda,k,post);

for k=1:4
ud(ud(:,:,:,k)==1)=k;
end

lab_n=reshape(permute(groundtruth,[2 1 3]),[],1);
lab_ud=reshape(permute(ud(:,:,:,1),[2 1 3]),[],1);
lab_ud=lab_ud-1;
error_segmentation = sum(lab_ud ~= lab_n) / numel(lab_n);
confusion_matrix = confusionmat(lab_ud,lab_n);
segmentation=ud(:,:,:,1);
end