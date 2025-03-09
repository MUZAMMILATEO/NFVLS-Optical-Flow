function [u v] = FSFOV_warp_2021(im1,im2,pyr_lev,warp_num,tempmat,wndcff,kernel_1,thetaOld,ep,lambdaFrac)
%THIS FUNCTION TAKES INPUT IMAGES AS IM1 AND IM2, 
%THE NUMBER OF PYRAMID LEVELS PYR_LEV, 
%THE NUMBER OF WARPING TIMES IN EACH PYRAMID LEVEL GIVEN BY WARP_NUM.

if size(im1,3)>1
   im1 = rgb2gray(im1);
   im2 = rgb2gray(im2);
end

if ~isfloat(im1)
   im1 = double(im1);
   im2 = double(im2);
end

[m n] = size(im1);

for lev = pyr_lev:-1:1

   tmp_im1 = imresize(im1,[ceil(0.5^(lev-1)*m) ceil(0.5^(lev-1)*n)]);
   tmp_im2 = imresize(im2,[ceil(0.5^(lev-1)*m) ceil(0.5^(lev-1)*n)]);
   [new_m new_n] = size(tmp_im1);
   if lev == pyr_lev
      u = zeros(new_m,new_n);
      v = zeros(new_m,new_n);
   end
   u = imresize(u,[new_m new_n]);
   v = imresize(v,[new_m new_n]);
   tmp_im1 = imwarp(tmp_im1,-u,-v,'true');
   for warping_number = 1:1:warp_num
      [fx,fy,ft] = computeDerivatives(tmp_im1,tmp_im2);
      fxx = fx.*fx;
      fyy = fy.*fy;
      fxy = fx.*fy;
      fxt = fx.*ft;
      fyt = fy.*ft;
      
      [delta_u,delta_v] = FSFOV_2021(fx,fy,ft,fxx,fyy,fxy,fxt,fyt,thetaOld,tempmat,wndcff,ep,kernel_1,lev,lambdaFrac);
      
      u = u + delta_u;
      v = v + delta_v;
      delta_u = medfilt2(delta_u,[5 5]);
      delta_v = medfilt2(delta_v,[5 5]);
      u = medfilt2(u,[9 9]);
      v = medfilt2(v,[9 9]);
      tmp_im1 = imwarp(tmp_im1,-delta_u,-delta_v,'true');
     end
end
figure, imshow(tmp_im1)
end