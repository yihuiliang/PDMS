function MaskAct = TrimapExpansion(I, Trimap, MaxIterExp)
% Trimap and MaskAct ----------
ITrimap = Trimap;
RmaskB = (ITrimap < 50);  RmaskF = (ITrimap > 200);
% LabelExpansion -------------------------------------------------
RMaskFExp=RmaskF ; RMaskBExp=RmaskB ;
ExpThr_U=9/256 ;
ExpThr_D=1/256 ;
ExpThrDist = ExpThr_U-ExpThr_D ;
%MaxIterExp=18 ;

%'Trimap expansion ...'
for i=1 : MaxIterExp
    ExpDist =i ;
    ExpThr =ExpThr_U - i* ExpThrDist / MaxIterExp ;
    [RMaskFExp ,RMaskBExp] = LabelExpansion (I, RMaskFExp, RMaskBExp , ExpDist, ExpThr)  ;
    
end

RmaskF =  RMaskFExp ;
RmaskB =  RMaskBExp ;

MaskAct = zeros(size(I, 1), size(I, 2)) ; % 0= inactive , 1= definit foreground ,% 5= definit background , 3= potential pixels
MaskAct = RmaskB * 5 + RmaskF ;   % get the marked F&B pixels
MaskAct(MaskAct==0)=3 ;

% 把mask调回人看得懂的数据
MaskAct(MaskAct == 1) = 255;
MaskAct(MaskAct == 5) = 0;
MaskAct(MaskAct == 3) = 128;
MaskAct = uint8(MaskAct);
end

