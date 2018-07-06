function y = L3DOneOutput_FGP(X)

% this function is the transpose of the L_FGP function
% it is related to the paper of Amir Beck and Teboulle

p  =  X - circshift(X,[-1,0,0]);
p  =  p(1:(size(X,1)-1),:,:)   ;
p(size(X,1),:,:) = 0;

q  =  X - circshift(X,[0,-1,0]);
q  =  q(:,1:(size(X,2)-1),:)   ;
q(:,size(X,2),:)=0;

r  =  X - circshift(X,[0,0,-1]);
r  =  r(:,:,1:(size(X,3)-1))   ;
r(:,:,size(X,3))=0;


y = [p(:),q(:),r(:)] ;