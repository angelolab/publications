function [yd,xd]=planeCoreg(CRchannel,massDS,countsAllSFilt,plotName)
% Calculate the x and y offset for each image
% CRchannel is the channel to be used for coregistration
% massDS is a dataset with the mass info
% countsAllSFilt is a cell array with length as the number of depths and in
% each cell there is a ImageXxImageYxMarkerNum matrix

yd=zeros(length(countsAllSFilt),1);
xd=zeros(length(countsAllSFilt),1);

[tf,channelLoc] = ismember(CRchannel,massDS.Label);

%Read in first plane as reference for all coregistrations
target=countsAllSFilt{1}(:,:,channelLoc);
Imagesize=size(target);

%Relative offset allowed, roff
%Maxiumium amount of offset in pixels that is allowed, s
roff=0.015;
sy=round(roff*Imagesize(1));
sx=round(roff*Imagesize(2))

%Calculate offsets
for i=1:length(countsAllSFilt)
    template=countsAllSFilt{i}(:,:,channelLoc);
    if ~(sum(sum(template))==0) & ~(sum(sum(target))==0)
        g=normxcorr2(template,target);
        gtrim=g(Imagesize(1)-sy:Imagesize(1)+sy,Imagesize(2)-sx:Imagesize(2)+sx);
    %   [yd(i),xd(i)]=find(g==max(g(:)));
        [yd(i),xd(i)]=find(gtrim==max(gtrim(:)),1,'First');
        yd(i)=yd(i)+Imagesize(1)-sy-1;
        xd(i)=xd(i)+Imagesize(2)-sx-1;
        disp(['Calculating coregistration offset for Depth ', num2str(i),'...','offset is xd=',num2str(xd(i)),', yd=',num2str(yd(i))])
    else
        disp(['No bits in image to calculate coregistration for ', CRchannel,' plane ', num2str(i)]);
        xd(i)=NaN;
        yd(i)=NaN;
    end
end

% plot
f1=figure;
subplot(2,1,1);
plot(xd-Imagesize(1),'*-');
xlabel('Depth');
ylabel('Offset');
title('Coregistration offset in X');

subplot(2,1,2);
plot(yd-Imagesize(2),'*-');
xlabel('Depth');
ylabel('Offset');
title('Coregistration offset in Y');
saveas(f1,plotName);
