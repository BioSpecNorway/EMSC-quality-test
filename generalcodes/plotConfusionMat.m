function plotConfusionMat(M1,cal_val)
% copied and edited from 
%https://stackoverflow.com/questions/3942892/how-do-i-visualize-a-matrix-with-colors-and-values-displayed

if ~iscell(M1)
    M{1} = M1;
else M = M1;
end

n = length(M);

for i = 1:n
sm = sum(M{i}.statistics.confusionMat,2);
Ngr = size(sm,1);
matrix = bsxfun(@rdivide,M{i}.statistics.confusionMat,sm);
figure
imagesc(matrix);
colormap(flipud(gray)); 
textStrings = num2str(matrix(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

[x,y] = meshgrid(1:Ngr);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(matrix(:) > midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gcf,'Color',[1 1 1]);
set(gca,'XTick',1:Ngr,...                         %# Change the axes tick marks
        'XTickLabel','',...
        'XAxisLocation', 'bottom',...%#   and tick labels
        'YTick',1:Ngr,...
        'YTickLabel',[M{i}.statistics.groupNames,...
        repmat('(',Ngr,1),num2str(sm),repmat(')',Ngr,1)],...
        'TickLength',[0 0]);
    ax = axis;
tx1 = text(1:Ngr,ax(3)*ones(1,Ngr),M{i}.statistics.groupNames);
set(tx1,'HorizontalAlignment','right','VerticalAlignment','top', ...
      'Rotation',-90,'FontSize',12,'FontWeight','bold');
    
if strcmp(cal_val,'cal')
    tx2 = text((ax(1)+ax(2))/2,ax(4),...
        ['Predicted class, Acc_{CV}=',num2str(100*M{i}.statistics.accuracy,'%0.1f'),'%']);
else
    tx2 = text((ax(1)+ax(2))/2,ax(4),...
         ['Predicted class, Acc_{Val}=',num2str(100*M{i}.statistics.accuracy,'%0.1f'),'%']);
end
set(tx2,'HorizontalAlignment','center','VerticalAlignment','top', ...
      'FontSize',12,'FontWeight','bold');

ty = text(ax(2),(ax(3)+ax(4))/2,'True class');
set(ty,'HorizontalAlignment','center','VerticalAlignment','top', ...
      'Rotation',90,'FontSize',12,'FontWeight','bold');
set(gca,'FontSize',12,'FontWeight','bold')
end

end
    
