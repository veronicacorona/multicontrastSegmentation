% Author: Veronica Corona, vc324@cam.ac.uk
% get labels and plot when testing a dataset

function [labels, error_classification]=getlabels_test(test,classifier, groundtruth)
% Compute local classification on test set using the specific classifier
labels=test*classifier*labeld;

d1=size(groundtruth,1);
d2=size(groundtruth,2);
d3=size(groundtruth,3);

% Plot labels
figure('Position',[10 10 588 549])
for s=1:d3
subplot(4,6,s)
imagesc(reshape(labels(((s-1)*d1*d2+1):(s*d1*d2))',d2,d1)');
axis off
colormap jet
end
% Misclassification error
lab_n=reshape(permute(groundtruth,[2 1 3]),[],1);
error_classification = sum(labels ~= lab_n) / numel(lab_n);
end

