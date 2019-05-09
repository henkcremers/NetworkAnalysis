function [dataRes] = nwa_reshape(data,method)
% ==================================================
% [dataRes] = nwa_reshape(data,method)
% ==================================================
%
% EXAMPLES
% 1) Switch between matrices and vectors
% mat = [1 2 3; 4 5 6; 7 8 9]
% dataRes = nwa_reshape(mat,'mat2vec') %returns the lower triangle values
% mat2 = nwa_reshape(dataRes,'vec2mat')
% 
% 2) Switch between 3D and 2D NWA structures
% mat(:,:,1) = [1 2 3; 4 5 6; 7 8 9]
% mat(:,:,2) = [10 11 12; 13 14 15; 16 17 18]
% mat(:,:,3) = [19 20 21; 22 23 24; 25 26 27]
% mat(:,:,4) = [28 29 30; 31 32 33; 34 35 36]
% [dataRes] = nwa_reshape(mat,'3d2d') % return a nn(3) x (n*(n-1))/2 matrix
% mat2 = nwa_reshape(dataRes,'2d3d')  % puts it back
% =================================================
nn = size(data);

%if ~isequal(nn(1),nn(2))

switch method
    
    case 'mat2vec'     
        
        L = tril(ones(nn(1),nn(1)),-1);
        dataRes = data(logical(L));
        
    case 'vec2mat'           
  
         dataRes = squareform(data);      
         
    case '3d2d'    

        dataRes = reshape(data,nn(1)^2,nn(3)); dataRes  = dataRes';
        L = tril(ones(nn(1),nn(1)),-1);
        lvec = reshape(L,1,nn(1)*nn(1)); 
        dataRes = dataRes(:,logical(lvec));
        
    case '2d3d'
        
        for j = 1:nn(1);           
            dataRes(:,:,j) = squareform(data(j,:));
        end
     
    case 'flipdiag'    
        dataRes   = rot90(fliplr(data),1);
        
end

return
