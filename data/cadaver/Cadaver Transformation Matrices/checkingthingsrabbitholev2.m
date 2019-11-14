% Checking slicer
% first load lps transform and reshape
clear all; clc;

rasnode = [0,0,1,1;0,1,0,2;-1,0,0,3;0,0,0,1];
lps2ras = diag([-1,-1,1,1]);
lpsnode = lps2ras*rasnode*lps2ras;
lpsfinal = inv(lpsnode);

load('C:\Users\riojaske\Desktop\TestNode.mat');
a = AffineTransform_double_3_3;
a = reshape(a,3,4)';
a = [a(1:3,1:3),a(4,:)'];
savedT_lps = [a;[0,0,0,1]]

% Initialize
itkoffset = zeros(3,1);
itkR = zeros(3,3);
vtkDim = 4;

for ii=1:(vtkDim-1)
    for jj=1:(vtkDim-1)
        itkR(ii,jj) = lpsnode(ii,jj);
    end
    itkoffset(ii) = lpsnode(ii,vtkDim);
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

% First Compute Offset lps
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