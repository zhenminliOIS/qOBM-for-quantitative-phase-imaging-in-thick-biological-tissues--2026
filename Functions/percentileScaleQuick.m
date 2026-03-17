function [cMin,cMax]=percentileScaleQuick(A,pMin,pMax,varargin)
%percentileScale(A,pMin,pMax,nBins) 
%Patrick Ledwig 1/18/19
%get max and min values for an upper and lower percentile of matrix for use
%when scaling an image.
%INPUTS:
%A - matrix
%pMin,pMax - percentiles (0-1) for cutoff
%nBins (optional) - number of bins for the histogram
%OUTPUTS:
%cMin,cMax - the ordinants of CDF where the respective percentiles fall
sA=sort(A(:));
iA=[1:numel(A(:))]./numel(A(:));
cMin=sA(find(iA>=pMin,1,'first'));
cMax=sA(find(iA>=pMax,1,'first'));
end