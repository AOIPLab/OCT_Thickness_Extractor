function [thickness] = funCalcThickness(segmentation, vertical_scale,lateral_scale, seed)
% This function is taken from the original calcThickness script. This
% converts the data from docTrap into scaled data and calculated the
% thicknesses of each layer.


% would need to adjust this to include any options for segmentations - JG
nPixels = length(segmentation(1,:));

tTRT=(segmentation(3,:)-segmentation(1,:)).*vertical_scale;
tIR=(segmentation(2,:)-segmentation(1,:)).*vertical_scale;
tOR=(segmentation(3,:)-segmentation(2,:)).*vertical_scale;
tChor=(segmentation(4,:)-segmentation(3,:)).*vertical_scale;
umThickness=[tTRT;tIR;tOR;tChor];

Xum=(1:nPixels).*lateral_scale;

cX=Xum-Xum(seed);

thickness(:,:,1)=[cX(1,:);umThickness(:,:,1)];



end