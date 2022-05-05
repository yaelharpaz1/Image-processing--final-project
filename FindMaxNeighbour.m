function [i_tag,j_tag] = FindMaxNeighbour(u,i,j)
%%
[m,n] = size(u);
  
if i~=1 && i~=m && j~=1 && j~=n
    A = [u(i,j) , u(i-1,j) , u(i,j-1) , u(i+1,j) , u(i,j+1)];
    index = [i,j, i-1,j, i,j-1, i+1,j, i,j+1];
    [~,I] = max(A);
    i_tag = index(2*I-1);
    j_tag = index(2*I);
else            
    if i==1 && j==1
        A = [u(i,j) , u(i+1,j) , u(i,j+1)];
        index = [i,j, i+1,j, i,j+1];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
    if i==1 && j==n
        A = [u(i,j) , u(i+1,j) , u(i,j-1)];
        index = [i,j, i+1,j, i,j-1];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
    if i==m && j==1
        A = [u(i,j) , u(i-1,j) , u(i,j+1)];
        index = [i,j, i-1,j, i,j+1];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
    if i==m && j==n
        A = [u(i,j) , u(i-1,j) , u(i,j-1)];
        index = [i,j, i-1,j, i,j-1];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
    if i==1 && j~=1 && j~=n
        A = [u(i,j) , u(i,j-1) , u(i+1,j) , u(i,j+1)];
        index = [i,j, i,j-1, i+1,j, i,j+1];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
    if i==m && j~=1 && j~=n
        A = [u(i,j) , u(i-1,j) , u(i,j-1) , u(i,j+1)];
        index = [i,j, i-1,j, i,j-1, i,j+1];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
    if j==1 && i~=1 && i~=m
        A = [u(i,j) , u(i-1,j) , u(i+1,j) , u(i,j+1)];
        index = [i,j, i-1,j, i+1,j, i,j+1];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
    if j==n && i~=1 && i~=m
        A = [u(i,j) , u(i-1,j) , u(i,j-1) , u(i+1,j)];
        index = [i,j, i-1,j, i,j-1, i+1,j];
        [~,I] = max(A);
        i_tag = index(2*I-1);
        j_tag = index(2*I);
    end
end

    
        
