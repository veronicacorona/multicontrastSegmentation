% Author: Veronica Corona, vc324@cam.ac.uk


% Compute labels, posterior probabilities and classification errors in test dataset
function [labels, error_classification, post, error_segmentation, confusion_matrix]=test_dataset(dataset, classifier, groundtruth,lambda, k)


%[test_labels, error_classification]=getlabels_test(test,classifier, groundtruth)
[labels, error_classification]=getlabels_test(dataset,classifier, groundtruth);
% Compute posterior maps
post=getposteriors_test(dataset,labels,classifier,groundtruth);

% Convex segmentation
%lambda=1;
error_segmentation = zeros (1,size(k,2));
%k=[0, 0.1, 0.2, 0.3, 0.4, 0.5];
for i=1:size(k,2)
[ud,u,error_segmentation(i),confusion_matrix(:,:,i)]=convex_segmentation(lambda,groundtruth,k(i), post);
end
error_segmentation=error_segmentation';
end
