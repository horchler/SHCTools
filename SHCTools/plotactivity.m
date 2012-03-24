function plotactivity(p,a,t,opts)
%PLOTACTIVITY  Plot simulated SHC network activity.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 4-8-10
%   Revision: 1.0, 3-24-12


figure
axis ij
n=size(a,1);
[y yt]=gety(p);
yy=cumsum(y);
yvec=[yy-y yy(end)]+0.5;
yaxis=yy(end)+0.5;
ytick=yy+0.5*(1-y);
xtick=0;
spacerows=0;
if nargin>3
    if ~iscell(opts)
        opts={opts};
    end
    for i=1:length(opts)
        switch opts{i}
            case 'noscaling'
                yvec=(0:n)+0.5;
                yaxis=n+0.5;
                ytick=1:n;
            case 'nodenumbers'
                yt=1:n;
            case 'timescale'
                xtick=1;
                set(gca,'XTickLabel',[])
            case 'spacerows'
                spacerows=1;
            case 'grayscale'
                colormap(1-gray)
        end
    end
end

if spacerows
    hold on
    %sf=[0 0.1;0.1 0.2;0.2 0.3;0.3 0.4;0 0.1;0.1 0.2;0.2 0.3;0.3 0.4];
    for i=1:n
        pcolor(t,yvec(i:i+1)-[0 0.01*(yaxis-0.5)],a([i i],:))
        %pcolor(t,yvec(i:i+1)-[0 0.01*(yaxis-0.5)]-sf(i,:),a([i i],:))
    end
    ytick=ytick-0.005*(yaxis-0.5);%-sf(:,2)'+0.05;
    yaxis=yaxis-0.01*(yaxis-0.5);%-0.4;
else
    pcolor(t,yvec,a([1:n n],:))
end
if xtick
    yaxis=1.025*yaxis;
    h=line([0 0.1*t(end)],[0 0]+0.999*yaxis);
    set(h,'Color',[0 0 0],'LineWidth',2,'LineSmoothing','on')
    text(0,1.005*yaxis,[num2str(0.1*t(end)) ' sec.'],'VerticalAlignment','top')
else
    xlabel('Time')
end
axis([t(1) t(end) 0.5 yaxis])
set(gca,'YTick',ytick,'YTickLabel',yt,'TickLength',[0 0])
set(gcf,'Color',[1 1 1])
box off
shading flat
ylabel('State ID')

function [y yt]=gety(p)
m=length(p);
y=zeros(1,m);
yt=cell(1,m);
j=1;
for i=1:m
    if iscell(p{i,i})
        [z zt]=gety(p{i,i});
        lz=length(z);
        y(j:j+lz-1)=z/length(p{i,i});
        for k=1:lz
            yt{j+k-1}=[num2str(i) ',' zt{k}];
        end
        j=j+lz;
    else
        y(j)=1;
        yt{j}=num2str(i);
        j=j+1;
    end
end