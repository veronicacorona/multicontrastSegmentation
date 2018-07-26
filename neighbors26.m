% Author: Veronica Corona, vc324@cam.ac.uk


function [m, sd]=neighbors26(matrix3d)
M1=padarray(matrix3d,[1 1 1]);
count=1;
for k=2:(size(M1,3)-1)
    for j=2:(size(M1,1)-1)
        for i=2:(size(M1,2)-1)
            w=M1(j-1:j+1,i-1:i+1,k-1:k+1);
            m(count)=mean(mean(mean(w)));
            for l=1:size(w,3)
            w1(l,:)=reshape(w(:,:,l)',1,9);
            end
            w2=reshape(w1',1,27);
            sd(count)=std(w2);
            count=count+1;
        end
    end
end
end
