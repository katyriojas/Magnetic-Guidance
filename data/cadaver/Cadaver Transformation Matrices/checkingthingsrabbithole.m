% Checking slicer
% first load lps transform and reshape
clear all; clc;

load('C:\Users\riojaske\Desktop\T_tracker_magsensor.mat');
a = AffineTransform_double_3_3;
a = reshape(a,3,4);
savedT_lps = [a;[0,0,0,1]]


lps2ras = diag([-1,-1,1,1]);

% Similarity Transform to express transformation in RAS slicer space
savedT_ras = lps2ras*savedT_lps*inv(lps2ras);

% Still need to handle the flip in order of rotation and translation
savedT_ras(1:3,4) = -savedT_ras(1:3,1:3)*savedT_ras(1:3,4);

% first itk converts to lps
savedT_lps = lps2ras*savedT_ras*lps2ras;
savedT_lpsInv = inv(savedT_lps);

% Initialize
itkoffset = zeros(3,1);
itkR = zeros(4,4);
vtkDim = 4;
%vtkDim is actually 4 (defined as 3 in c++ zeros base)

for ii=1:(vtkDim-1)
    for jj=1:(vtkDim-1)
        itkR(ii,jj) = savedT_lpsInv(ii,jj);
    end
    itkoffset(ii) = savedT_lpsInv(ii,vtkDim);
end
 
% These are what are fed into the affine transform:
% itkR 
% itkoffset

% Initialize things, only thing that should change is the offset

% I think this should actually be three because this is setting size of the
% array
m_cent = zeros(vtkDim-1,1);
m_Trans = zeros(vtkDim-1,1);
offset = zeros(vtkDim-1,1);
translation = zeros(vtkDim-1,1);

% First Compute Offsetlps
for ii = 1:vtkDim-1
     
    offset(ii) = m_Trans(ii) + m_cent(ii);
    
    for jj = 1:vtkDim-1
        
        offset(ii) = m_cent(ii) - itkR(ii,jj)*m_cent(jj)
        
    end
    
end

m_Offset1 = offset; % this should be zero because just rotational part updated

m_Offset = itkoffset; % set the offset with the translational comp of rigid transform

% Now Compute Translation
for ii = 1:(vtkDim-1)
    
    translation(ii) = m_Offset(ii) - m_cent(ii);
    
    for jj=1:(vtkDim-1)
        
        translation(ii) = translation(ii) + itkR(ii,jj)*m_cent(jj);
        
    end
    
end

itkR
m_Translation = translation