% Author: Veronica Corona, vc324@cam.ac.uk

function [n1, n2, n3, n4, n5, n6]=find_six_neighbors(matrix3d)

M=padarray(matrix3d,[1 1 1]);
n1=zeros(1,numel(matrix3d));
n2=zeros(1,numel(matrix3d));
n3=zeros(1,numel(matrix3d));
n4=zeros(1,numel(matrix3d));
n5=zeros(1,numel(matrix3d));
n6=zeros(1,numel(matrix3d));

count=1;
for k=2:(size(M,3)-1)
    for j=2:(size(M,1)-1)
        for i=2:(size(M,2)-1)
            
            n1(count)=M(j,i-1,k);
            n2(count)=M(j,i+1,k);
            n3(count)=M(j-1,i,k);
            n4(count)=M(j+1,i,k);
            n5(count)=M(j,i,k-1);
            n6(count)=M(j,i,k+1);
            count=count+1;
        end
    end
end
end
