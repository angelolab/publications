function obj=applyOffset(xd,yd,depth,Imagesize,FilteredChannel)
%  specify start and stop indices for rows and columns
% [i,j,xd(j),yd(j)]
CurrentFrame=uint16(zeros(Imagesize));
xoff=Imagesize(2)-xd(depth);
         yoff=Imagesize(1)-yd(depth);
%          
         if yoff<0;
             ai=abs(yoff)+1;
             bi=Imagesize(1);
             ci=1;
             di=Imagesize(1)-abs(yoff);
         end
         if yoff>0;
             ai=1;
             bi=Imagesize(1)-abs(yoff);
             ci=abs(yoff)+1;
             di=Imagesize(1);
         end
         if yoff==0;
             ai=1;
             bi=Imagesize(1);
             ci=1;
             di=Imagesize(1);
         end
         if xoff<0;
             ei=abs(xoff)+1;
             fi=Imagesize(2);
             gi=1;
             hi=Imagesize(2)-abs(xoff);
         end
         if xoff>0;
             ei=1;
             fi=Imagesize(2)-abs(xoff);
             gi=abs(xoff)+1;
             hi=Imagesize(2);
         end
         if xoff==0;
             ei=1;
             fi=Imagesize(2);
             gi=1;
             hi=Imagesize(2);
         end

         CurrentFrame(ai:bi,ei:fi)=FilteredChannel(ci:di,gi:hi);
         obj=CurrentFrame;
end
        
