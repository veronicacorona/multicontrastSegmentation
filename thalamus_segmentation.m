%% A multi-contrast MRI segmentation appraoch to thalamus segmentation


% Author Information: 
% Veronica Corona
% Department of Applied Mathematics and Theoretical Physics
% University of Cambridge, UK

% Contact: vc324@cam.ac.uk

% We propose a method to automatic segment the posterior, lateral and medial
% nucleus in the thalamus. The algorithm constists of two stages: first, we perform
% a supervised learning approach to obtain classification of each voxel;
% second, we refine the segmentation solving a convex optimization problem.
% In particulare we solve: u = argmin_u - log P(u|I) + \lambda TV(u). 

% This code is implemented in MATLAB to perform multi-contrast segmentation as 
% described in the preprint of the article at http://arxiv.org/. Here, we present 
% the application of thalamic nuclei, but this code can be used for any other multi-contrast
% dataset for which you have manual segmentations. The feature space can also 
% be tuned differently for different applications. 
% If you use this code in your research, please cite this work. 

warning off
prmemory(200000000000000000000)
prwaitbar off
%% Load dataset43 and its manual segmentation
% We select a region of interest containing the thalamus from dataset43. We
% use this dataset as training set. We use the manual segmentation of the major nuclei
% as ground truth in the supervised learning stage.

% Load training set dataset 43
QSM_mean=load_nii('all_QSM_mean.nii'); 
MPRAGEonT2s=load_nii('MPRAGEonT2s.nii');
T2s_template=load_nii('T2star_magnitude_template.nii'); 
% load manual segmentation for dataset 43
nuclei43=load_nii('Thalamic_Nuclei.nii'); 

% Select region of interest
t2_left=double(rot90(select_left_thalamus(T2s_template.img)));
m_left=double(rot90(select_left_thalamus(MPRAGEonT2s.img)));
q_left=double(rot90(select_left_thalamus(QSM_mean.img)));

n_left=double(rot90(select_left_thalamus(nuclei43.img)));
% Load manual segmentation as labels for training
mask43=load_nii('thalamus43.nii.gz');
mask=double(rot90(select_left_thalamus(mask43.img)));
% Create a vector of labels that will be used in our voxel-based dataset
labels43=reshape(permute(mask,[2 1 3]),[],1);

%% Create a voxel-based dataset.
% Each voxel represents an instance in the dataset. We compute 27 features in total: 3 constrasts intensity values, mean
% and std in the 26-neighboors, intensity values of 6 closest neighbors.

% Create dataset using the 3 contrasts and 27.
num_feat = 27; % It is possible to change the set of features for a different problem/dataset. We implemented a few alternatives
dataset=create_dataset(t2_left,m_left,q_left,num_feat);
% assign a label to each voxel
data=prdataset(dataset,labels43);

% Set equal prior probability for each class in the dataset
prior(1)=size(find(labels43==0),1)/numel(q_left);
prior(2)=size(find(labels43==1),1)/numel(q_left);
prior(3)=size(find(labels43==2),1)/numel(q_left);
prior(4)=size(find(labels43==3),1)/numel(q_left);
data=setprior(data,prior);

% Scale features by variance
data=data*scalem(data,'variance');
%% Training stage
% We are using k Nearest Neighbors (k=3) and Parzen Classifier
% nn3, parzen_cl are objects that can be used to classify a new test set:
% labels=test*classifier*labeld;
% to otain posterior probabilities: posterior=testset*classifier*classc;

nn3=knnc(data,3);
[parzen_cl, h]=parzenc(data); % h is the optimal width of the window density estimator which is computed automatically
%% Test templates from dataset116
% Load dataset116
mmprage=load_nii('avg_all-age_mMPRAGE.nii.gz');
avg_qsm=load_nii('avg_all-age_ref01_QSM01.nii.gz');
n4_magn=load_nii('n4_avg_all-age_n4_magn.nii.gz');
% Select the region of interest for the three contrasts
MPRAGE=double(rot90(select_left_thalamus_small(mmprage.img)));
T2s=double(rot90(select_left_thalamus_small(n4_magn.img)));
QSM116=double(rot90(select_left_thalamus_small(avg_qsm.img)));
% Load manual segmentation that will be used to validate the result of the automatic
% segmentation
nuclei = load_nii('thalamus_all.nii.gz');
n116=double(rot90(select_left_thalamus_small(nuclei.img)));
% this adjusts values in the mask to give class labels as follows:
% 0 = background
% 1 = lateral nucleus
% 2 = medial nucleus
% 3 = posterior nucleus
n116(n116==4)=3;
n116(n116==5)=1;

% Create dataset with 27 features as we did for the training set
dataset1=create_dataset(T2s,MPRAGE,QSM116,27);
dataset116=prdataset(dataset1);
dataset116=dataset116*scalem(dataset116,'variance');
%% Compute labels, posterior probabilities and classification errors in dataset116
% Weight for prior knowledge of nuclei locations
K=linspace(0,1,11);

% Preallocation
error_segmentation_nn3 = zeros(size(K,2),1);
error_segmentation_nn5 = zeros(size(K,2),1);
error_segmentation_p = zeros(size(K,2),1);
confusion_matrix_nn3 = zeros(4,4,size(K,2));
confusion_matrix_nn5 = zeros(4,4,size(K,2));
confusion_matrix_p = zeros(4,4,size(K,2));

%% 3NN classifier
% Compute labels, posterior probabilities and classification errors in dataset116

% Compute labels using 3NN of dataset 116 and compute the local
% classification error as number of misclassifed voxels.
[labels116_nn3, error_classification(1)] = getlabels_test(dataset116,nn3, n116);
% Compute posterior maps
post_nn3 = getposteriors_test(dataset116,labels116_nn3,nn3);
% Convex segmentation to refine the segmentation: we produce a 3D segmentation
% ud, misclassification error and confusion matrix for the final classification.
lambda=1; % regularisation parameter (TO DO: show that it is consistent with the choice of classifier
for i=1:size(K,2)
    [ud_nn3(:,:,:,i),error_segmentation_nn3(i),confusion_matrix_nn3(:,:,i)]=convex_segmentation(lambda,n116,K(i), post_nn3);
end

%% Parzen classifier
% Compute labels, posterior probabilities and classification errors in dataset116
[labels116_p, error_classification(2)] = getlabels_test(dataset116,parzen_cl, n116);
% Compute posterior maps
post_p = getposteriors_test(dataset116,labels116_p,parzen_cl);
% Convex segmentation to refine the segmentation: we produce a segmentation
% ud, misclassification error and confusion matrix for the final classification.
lambda=5.5;
for i=1:size(K,2)
    [ud_parz(:,:,:,i),error_segmentation_p(i),confusion_matrix_p(:,:,i)]=convex_segmentation(lambda,n116,K(i), post_p);
end
%% Classifiers comparison
% Local classification error vs convex segmentation error (for k=0)
error_segmentation=[error_segmentation_nn3(1) error_segmentation_p(1)];
EE=[error_classification' error_segmentation'];
figure
bar(EE);
Labels = {'3NN', 'Parzen'};
set(gca, 'XTick', 1:4, 'XTickLabel', Labels);
ylabel('Error')
legend('Classification error','Convex segmentation error')
%% Load Morel Atlas
load_morelatlas;
% Compute error of the morel atlas segmentation on the 116 dataset
error_morel = size(find(atlas ~= n),1 )/ numel(n); % n is the groundtruth for dataset 116

%% Single subjects (patients 10, 39, 80 with manual segmentations) varying k
patient=10; % try 10,39, 80
[T2s,MPRAGE,QSM,n]=load_single(patient);
lab_n10=reshape(permute(n,[2 1 3]),[],1);

dataset2=create_dataset_new(T2s,MPRAGE,QSM,27);
dataset22=prdataset(dataset2);
datanew1=dataset22*scalem(dataset22,'variance');

K=linspace(0,1,11);
lambda=1;
classifier = nn3;
[~, error_classification_nn3, ~, error_segmentation_nn3,~]=test_dataset(datanew1, classifier, n,lambda, K);

lambda=5;
classifier = parzen_cl;
[~, error_classification_parz, ~, error_segmentation_parz,~]=test_dataset(datanew1, classifier, n,lambda, K);
