function[p]= plotPotts(so,s1,s2,s)

subplot(2,2,1)
imagesc(so);axis equal;axis off;
colormap('default');

subplot(2,2,2)
imagesc(s1),axis equal;axis off;

subplot(2,2,3)
imagesc(s2);axis equal;axis off;    

subplot(2,2,4)
imagesc(s);axis equal;axis off;    

