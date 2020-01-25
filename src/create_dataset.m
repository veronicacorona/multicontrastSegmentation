% Author: Veronica Corona, vc324@cam.ac.uk

% Given the three contrasts, this function builds the dataset in a voxel
% based representation, containing num_feat features. We use 27 features,
% but other possibilities are implemented
function [dataset]=create_dataset(t2_left,m_left,q_left,num_feat)
% 3 contrast intensity values
dataset(:,1)=reshape(permute(t2_left,[2 1 3]),1,numel(t2_left))';
dataset(:,2)=reshape(permute(m_left,[2 1 3]),1,numel(t2_left))';
dataset(:,3)=reshape(permute(q_left,[2 1 3]),1,numel(t2_left))';

% 6 closest neighbours
[n1_t, n2_t, n3_t, n4_t, n5_t, n6_t]=find_six_neighbors(t2_left);
[n1_m, n2_m, n3_m, n4_m, n5_m, n6_m]=find_six_neighbors(m_left);
[n1_q, n2_q, n3_q, n4_q, n5_q, n6_q]=find_six_neighbors(q_left);
dataset(:,4)=n1_t;
dataset(:,5)=n2_t;
dataset(:,6)=n3_t;
dataset(:,7)=n4_t;
dataset(:,8)=n5_t;
dataset(:,9)=n6_t;
dataset(:,10)=n1_m;
dataset(:,11)=n2_m;
dataset(:,12)=n3_m;
dataset(:,13)=n4_m;
dataset(:,14)=n5_m;
dataset(:,15)=n6_m;
dataset(:,16)=n1_q;
dataset(:,17)=n2_q;
dataset(:,18)=n3_q;
dataset(:,19)=n4_q;
dataset(:,20)=n5_q;
dataset(:,21)=n6_q;
if num_feat==21
else
    if num_feat>21
        [m_t, sd_t]=neighbors26(t2_left);
        [m_m, sd_m]=neighbors26(m_left);
        [m_q, sd_q]=neighbors26(q_left);
        
        % mean and std in the 26-neighbourhood
        dataset(:,22)=m_t;
        dataset(:,23)=sd_t;
        dataset(:,24)=m_m;
        dataset(:,25)=sd_m;
        dataset(:,26)=m_q;
        dataset(:,27)=sd_q;
       
        if num_feat==27
            
        else
            % histogram values in 16 bins of QSM
            [numelem]=hist_feat(q_left);
            dataset(:,28:(28+15))=numelem;
            
            for i=1:size(q_left,3)
                [Gmag(:,:,i), ~]=imgradient(q_left(:,:,i));
            end
            % histogram values in 16 bins for the gradient magnitude of QSM
            [n_gr]=hist_feat(Gmag);
            dataset(:,44:(44+15))=n_gr;
            
            % figure; montage_roi(gLap)
            
            % GAUSSIAN SPACE
            % clear h
            %
            % h=fspecial('gaussian',[3 3],1);
            % for i=1:size(q_left,3)
            % G1(:,:,i)=imfilter(q_left(:,:,i),h);
            % end
            %
            % h2=fspecial('gaussian',[5 5],2);
            % for i=1:size(q_left,3)
            % G2(:,:,i)=imfilter(q_left(:,:,i),h2);
            % end
            %
            % h3=fspecial('gaussian',[9 9],4);
            % for i=1:size(q_left,3)
            % G3(:,:,i)=imfilter(q_left(:,:,i),h3);
            % end
            % %
            % dataset(:,60)=reshape(permute(G1,[2 1 3]),[],1);
            % dataset(:,61)=reshape(permute(G2,[2 1 3]),[],1);
            % dataset(:,62)=reshape(permute(G3,[2 1 3]),[],1);
            
            % LAPLACIAN
            h=fspecial('laplacian');
            for i=1:size(q_left,3)
                gLap(:,:,i)=imfilter(q_left(:,:,i),h);
            end
            dataset(:,60)=reshape(permute(gLap,[2 1 3]),[],1);
        end
        
    end
end


end

