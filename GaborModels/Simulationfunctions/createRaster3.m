function [] = createRaster3(cellStats, preStimTime, postStimTime, color)
% inputs: cellStats- a matrix where the elements of each row are spike
% times from a different trial
% preStimTime - time before 0 to display raster
% postStimTime - time after 0 to display raster
% color - optional input of lenght equal to numTrials with a character for
% color

%scavenged from LstStitchEnsRaster02 from EST at DalyLab summer 2008
% plots raster in current figure window
if nargin < 4
    color = 'b';
end
APTraces = cellStats;
linWid = 2;
cumLen = size(cellStats,1);
y = [0.001,1];
for i = 1:size(cellStats,1)
    sp = [];
    sp = cellStats(i,cellStats(i,:)>0);
    if ~isempty(sp)
        for j = 1:length(sp) 
            if size(color,1)==1
                plot(ones(size(y))*sp(j), i+y,'color',color)
            else
                plot(ones(size(y))*sp(j), i+y,'color',color(i,:))                
            end
            hold on
        end
    end
    
end
set(gca,'YDir','reverse');
set(gca,'TickDir','out') % draw the tick marks on the outside
set(gca,'YTick', []) % don't draw y-axis ticks
box off
axis([-preStimTime postStimTime 0 cumLen+1 ]);
