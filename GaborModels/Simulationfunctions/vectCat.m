function [newMat] = vectCat(mat1, mat2)

[row1 col1] = size(mat1);
[row2 col2] = size(mat2);

newMat = zeros(max(row1, row2), col1+col2);
newMat(1:row1,1:col1) = mat1;
newMat(1:row2,col1+1:col1+col2) = mat2;

