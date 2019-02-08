function [stackLayers, peakPos] = sumSlicesPeak(model,direction,atomAddup)
% to return the summed layers between atomic layer spacings
% input:
%       model is the reconstructed model with same xzy dimensions as N
%       dimension = 1,2,3 (indicating viewing along rows, columns, pages)
% output:
%       stackLayers is the 3D atomic layers, all in N*N*m format.
% 2017.05.27 Xuezeng.
% Input best fit manually


[dim1 dim2 dim3]= size(model);


if nargin == 2
    atomAddup = 1; % how many neighboring pixels for an atomic layer
end

fprintf('Adding %d nearby pixels \n',atomAddup);

switch direction
    case 1
        sum1D=double(squeeze(sum(sum(model,2),3)));
        if size(sum1D,2)==1
            sum1D=sum1D';
        end
        peakPos = find(sum1D>circshift(sum1D,[0 1]) & sum1D>circshift(sum1D,[0 -1]));
        peakPos(find(peakPos==1))=[];
%         peakPos(find(peakPos==2))=[];
        peakPos(find(peakPos==dim1))=[];
        peakPos(find(peakPos==dim1-1))=[];
        peakPos(find(peakPos==dim1-2))=[];
        figure();plot(1:1:dim1,sum1D,'b-',peakPos,sum1D(peakPos),'r.','MarkerSize',14);title('sumX');
        
        stackLayers = zeros(dim2,dim3);
        for i = peakPos
            stack=squeeze(sum(model(i-atomAddup:i+atomAddup,:,:),direction));
            stackLayers = cat(3,stackLayers,stack);
        end 
        stackLayers(:,:,1)=[];
    case 2
        sum1D=double(squeeze(sum(sum(model,1),3)));
        if size(sum1D,2)==1
            sum1D=sum1D';
        end
        peakPos = find(sum1D>circshift(sum1D,[0 1]) & sum1D>circshift(sum1D,[0 -1]));
        peakPos(find(peakPos==1))=[];
        peakPos(find(peakPos==2))=[];
        peakPos(find(peakPos==dim2))=[];
        peakPos(find(peakPos==dim2-1))=[];
        figure();plot(1:1:dim2,sum1D,'b-',peakPos,sum1D(peakPos),'r.','MarkerSize',14);title('sumY');
        
        stackLayers = zeros(dim1,dim3);
        for i = peakPos
            stack=squeeze(sum(model(:,i-atomAddup:i+atomAddup,:),direction));
            stackLayers = cat(3,stackLayers,stack);
        end 
        stackLayers(:,:,1)=[];
        
    case 3
        sum1D=double(squeeze(sum(sum(model,1),2)));
        if size(sum1D,2)==1
            sum1D=sum1D';
        end
        peakPos = find(sum1D>circshift(sum1D,[0 1]) & sum1D>circshift(sum1D,[0 -1]));
        peakPos(find(peakPos==1))=[];
        peakPos(find(peakPos==2))=[];
        peakPos(find(peakPos==dim3))=[];
        peakPos(find(peakPos==dim3-1))=[];
        figure();plot(1:1:dim3,sum1D,'b-',peakPos,sum1D(peakPos),'r.','MarkerSize',14);title('sumZ');
        
        stackLayers = zeros(dim1,dim2);
        for i = peakPos
            stack=squeeze(sum(model(:,:,i-atomAddup:i+atomAddup),direction));
            stackLayers = cat(3,stackLayers,stack);
        end 
        stackLayers(:,:,1)=[];
end
end