% Author: Veronica Corona, vc324@cam.ac.uk

% COMPUTE POSTERIOR MAPS FOR EACH CLASS 

function [post]=getposteriors_test(test,labels_test,classifier,gt)

d1=size(gt,1);
d2=size(gt,2);
d3=size(gt,3);


% calculate posterior probabilities
data2 = prdataset(test,labels_test);
p = data2*classifier*classc;
pp = +p+0.000001; % avoid zeros
clear post
post = zeros(36,25,24,4);
% create the 39x30x24x4 box to do 3D convex segmentation
for k=1:24
    post(:,:,k,1)=reshape(pp(((k-1)*d1*d2+1):k*900,1)',d2,d1)';
    post(:,:,k,2)=reshape(pp(((k-1)*d1*d2+1):k*900,2)',d2,d1)';
    post(:,:,k,3)=reshape(pp(((k-1)*d1*d2+1):k*900,3)',d2,d1)';
    post(:,:,k,4)=reshape(pp(((k-1)*d1*d2+1):k*900,4)',d2,d1)';
end
save ('post3d.mat','post','labels_test','pp');

k=12; % plot posterior map for one slice only
figure
s(1) = subplot(1,4,1); imagesc(post(:,:,k,1), [0 1]); axis off; axis image; title('Background')
s(2) = subplot(1,4,2);imagesc(post(:,:,k,2), [0 1]); axis off; axis image; title('Lateral')
s(3) = subplot(1,4,3);imagesc(post(:,:,k,3), [0 1]); axis off; axis image; title('Medial')
s(4) = subplot(1,4,4); imagesc(post(:,:,k,4), [0 1]); axis off; axis image; title('Posterior')
colormap jet 
% Get positions of all the subplot
posa = get(s,'position');
h    = colorbar;
% Reset ax(3) position to before colorbar
set(s(4),'position',posa{4})
% Set everything to units pixels (avoids dynamic reposition)
set([s h],'units','pix');
% Widen figure by a factor of 1.1 (tweak it for needs)
posf = get(gcf,'position');
set(gcf,'position',[posf(1:2) posf(3)*1.1 posf(4)])

end
