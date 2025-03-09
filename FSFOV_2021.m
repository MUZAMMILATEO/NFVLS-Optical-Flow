function [u,v] = FSFOV_2021(fx,fy,ft,fxx,fyy,fxy,fxt,fyt,thetaOld,tempmat,wndcff,ep,kernel_1,lev,lambdaFrac)
%THIS FUNCTION ESTIMATES OPTICAL FLOW BASED ON FIELD SEGMENTATION
%FRACTIONAL ORDER VARIATIONAL MODEL.

[mIm,nIm] = size(fx);
u = zeros(mIm,nIm);
v = u;
u_cap = v;
v_cap = u_cap;
img = imread('C:\Users\acer\Desktop\GitHub_Code\Expert_System\circle.jpg');
a = round(5^(0.5*(lev-1)));
img = imresize(img,[a,a]);
img = mat2gray(img);
img = (img==1);
img = double(img);

%%%%%%%%%%%%%%%%%%%%%%% INITIALIZING THE LEVEL SURFACES %%%%%%%%%%%%%%%%%%%
phi_u = zeros(mIm,nIm);
numLevX = floor(mIm/a);
numLevY = floor(nIm/a);%floor(numLevX*nIm/mIm);
numLevX0 = numLevX-1;
numLevY0 = numLevY-1;

for i = 0:1:numLevX0
    for j = 0:1:numLevY0
        phi_u((a*i+1):a*(i+1),(a*j+1):a*(j+1)) = img;
    end
end

remX = mod(mIm,a);
remY = mod(nIm,a);

for i = 0:1:numLevX0
    phi_u((a*i+1):a*(i+1),a*numLevY+1:end) = img(:,1:remY);
end

for j = 0:1:numLevY0
    phi_u(a*numLevX+1:end,(a*j+1):a*(j+1)) = img(1:remX,:);
end

phi_u(a*numLevX+1:end,a*numLevY+1:end) = img(1:remX,1:remY);

 phi_v = phi_u;
% figure,imshow(phi_u)
%%%%%%%%%%%%%%%%%%%%%%% INITIALIZED THE LEVEL SURFACES %%%%%%%%%%%%%%%%%%%%

for it = 1:1:100
    [u_cap, v_cap, thetaOld] = FSFOV_2021_first(fx,fy,ft,fxx,fyy,fxy,fxt,fyt,u,v,u_cap,v_cap,thetaOld);
    u_cap = conv2(u_cap,kernel_1,'same');
    v_cap = conv2(v_cap,kernel_1,'same');
    [u_plus, u_minus, v_plus, v_minus, ep] = FSFOV_2021_second(u_cap,v_cap,u,v,thetaOld,phi_u,phi_v,tempmat,wndcff,ep,lambdaFrac);
    u_plus = conv2(u_plus,kernel_1,'same');
    u_minus = conv2(u_minus,kernel_1,'same');
    v_plus = conv2(v_plus,kernel_1,'same');
    v_minus = conv2(v_minus,kernel_1,'same');
    [phi_u, phi_v] = FSFOV_2021_third(phi_u,phi_v,u_plus,u_minus,v_plus,v_minus,u_cap,v_cap,ep,tempmat,thetaOld,lambdaFrac,lev);
    %[phi_u, phi_v] = Lu_2021_third(phi_u,phi_v,u_plus,u_minus,v_plus,v_minus,u_cap,v_cap,ep,thetaOld,lambdaFrac,lev);
    [u,v] = FSFOV_2021_fourth(u_plus,u_minus,v_plus,v_minus,phi_u,phi_v,ep);
    u = conv2(u,kernel_1,'same');
    v = conv2(v,kernel_1,'same');

end

end

