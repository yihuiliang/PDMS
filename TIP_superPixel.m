function [] = TIP_superPixel(png_name)
% Pixel-level Discrete Multiobjective Sampling Demo
% Yihui Liang, IAIS, SCUT
% 2019.07.20
% This demo implements the PDMS-based (pixel-level discrete
% multiobjective sampling) matting algorithm. Notice that preprocessing and
% postprocessing are not applied in this demo.
% Reference: Huang H, Liang Y*, Yang X, et al. Pixel-level Discrete
% Multiobjective Sampling for Image Matting[J]. IEEE Transactions on Image
% Processing, 27(8):3739-3751

image_name = strsplit(png_name,'.');
% feature_name = strcat('img_feature/img_feature_',image_name{1,1},'.mat');
% s = load(feature_name);
% img_feature = s.img_feature;
% img_feature = reshape(img_feature,size(img_feature,1)*size(img_feature,2),100);

addpath('expansion');
beginTime = clock;
out_dir = 'result\';
% png_name = 'GT01.png';
%% Load image and trimap
raw_img = imread(png_name);
trimap = imread(['Trimap2\',png_name]);
trimap = TrimapExpansion(double(raw_img), trimap, 11);

cluster_num = 150;

%% Pixel-level Discrete Multiobjective Sampling
%feature
U_ind = find(trimap ==128);
F_ind = find(trimap ==255);
B_ind = find(trimap ==0);
img = reshape(raw_img,numel(trimap),3);
F_color = img(F_ind,:);
B_color = img(B_ind,:);

% F_feature = img_feature(F_ind,:);
% B_feature = img_feature(B_ind,:);

U_rgb = img(U_ind,:);
[U_y,U_x] = ind2sub(size(trimap),U_ind); U_s  = [U_y,U_x];
% opts = statset('MaxIter',500); 
% idx = kmeans(single([U_rgb,U_s]), cluster_num,'Options',opts);
% save([out_dir,'idx_',png_name(1:4),'.mat'],'idx');
% idx = load([out_dir,'idx_',png_name(1:4),'.mat']);idx = idx.idx;
[labels, numlabels] = slicmex(raw_img,500,20);  %numlabels is the same as number of superpixels

[y,x] = ind2sub(size(trimap),1:numel(trimap));
U_alpha = zeros(length(U_ind),1);

step = 1;
for c = 0:numlabels-1
   cc_U_ind = find(labels==c);
%    ss = find(U_ind==cc_U_ind);
   [~,uind] = ismember(cc_U_ind, U_ind, 'rows');
   ind_of_U_ind = nonzeros(uind);
   if ~isempty(ind_of_U_ind)
      t1 = clock;
      n = ind_of_U_ind(randperm(length(ind_of_U_ind),1));
      U_color = single(img(U_ind(n),:));
%       
%       U_feature = img_feature(U_ind(n),:);     
      %foreground color sampling
      spatial_cost_FU = pdist2([x(F_ind)',y(F_ind)'],[x(U_ind(n)),y(U_ind(n))]);
      color_cost_FU = pdist2(single(F_color),single(img(U_ind(n),:)));
      pareto_F_bw = FDMO(double([spatial_cost_FU,color_cost_FU]));
      
      spatial_cost_BU = pdist2([x(B_ind)',y(B_ind)'],[x(U_ind(n)),y(U_ind(n))]);
      color_cost_BU = pdist2(single(B_color),single(img(U_ind(n),:)));
      pareto_B_bw = FDMO(double([spatial_cost_BU,color_cost_BU]));
      %sample feature extration 
      F_pareto_color = single(F_color(pareto_F_bw,:)');
      B_pareto_color = single(B_color(pareto_B_bw,:)');
%       F_pareto_feature = F_feature(pareto_F_bw,:)';
%       B_pareto_feature = B_feature(pareto_B_bw,:)';
      F_pareto_dist = spatial_cost_FU(pareto_F_bw);
      B_pareto_dist = spatial_cost_BU(pareto_B_bw);
      %evaluation and alpha estimation
      [U_alpha(ind_of_U_ind)] = EvaluateSamples(pareto_F_bw,pareto_B_bw,...
          F_pareto_color,B_pareto_color,F_pareto_dist,B_pareto_dist,U_color);
      t2=clock;
      e = etime(t2,t1);
      e = fprintf('%s ÆÀ¹ÀµÚ%d/%d¸ö³¬ÏñËØ¿é£¨ÄÚº¬%d¸öÎ´ÖªÏñËØ£©ºÄÊ±%fÃë\n',png_name(1:4),step,numlabels,length(ind_of_U_ind),e);
      step=step+1;
   end
end

alpha = trimap;
alpha(U_ind) = U_alpha*255;
endTime = clock;
fid = fopen([out_dir,'time.txt'],'a');
fprintf(fid, '%s total time:%f sec\n',png_name(1:4),etime(endTime,beginTime));
fclose(fid);
imwrite(alpha, [out_dir,png_name]);
clear,close all;
end
%  Fast Discrete Multiobjective Optimization Function
%Input:an array of all feasible solutions of the MOP denoted as C
%Output:C,i(the first i-1 elements in array C is Pareto optimal.In a row as
%a combination of solutions) and pareto_front_bw (the binary vector that
%indicates the Pareto front)
% function[pareto_front_bw,C,i] = FDMO(C)
% i=1;
% j=size(C,1);
% solution_ind = 1:size(C,1);
% while(i<=j)
%     cmp=i+1;
%     while(cmp<=j)
%         if( all(C(i,:)<=C(cmp,:))&&any(C(i,:)<C(cmp,:)) )
%             solution_ind([j,cmp])=solution_ind([cmp,j]);
%             C([j,cmp],:)=C([cmp,j],:);
%             j=j-1;
%         elseif(all(C(cmp,:)<=C(i,:))&&any(C(cmp,:)<C(i,:)))
%             C([i,cmp,j],:)=C([cmp,j,i],:);
%             solution_ind([i,cmp,j])=solution_ind([cmp,j,i]);
%             j=j-1;
%             cmp=i+1;
%         else
%             cmp=cmp+1;
%         end
%     end
%     i=i+1;
% end
% pareto_front_bw = false(size(C,1),1);
% pareto_front_bw(solution_ind(1:i-1)) = 1;
% end