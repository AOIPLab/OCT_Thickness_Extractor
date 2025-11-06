function [thickness] = calcThickness(segmentation, vertical_scale,lateral_scale)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nPixels = length(segmentation(1,:));

tTRT=(segmentation(3,:)-segmentation(1,:)).*vertical_scale;
tIR=(segmentation(2,:)-segmentation(1,:)).*vertical_scale;
tOR=(segmentation(3,:)-segmentation(2,:)).*vertical_scale;
tChor=(segmentation(4,:)-segmentation(3,:)).*vertical_scale;
umThickness=[tTRT;tIR;tOR;tChor];

Xum=(1:nPixels).*lateral_scale;

cX=Xum-Xum(XONH);

thickness(:,:,1)=[cX(1,:);umThickness(:,:,1)];



end