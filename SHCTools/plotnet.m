function plotnet(p,opts)
%PLOTNET  Visualization of SHC network topology.
%
%

%   Andrew D. Horchler, adh9@case.edu, Created 4-2-10
%   Revision: 1.0, 3-24-12


grayscale=0;
selfexcitation=1;
nodenumbers=0;
s='';
if nargin>1
    if ~iscell(opts)
        opts={opts};
    end
    for i=1:length(opts)
        switch opts{i}
            case 'noselfexcitation'
                selfexcitation=0;
            case 'grayscale'
                grayscale=1;
            case 'nosubnets'
                [p s]=buildrho(p);
                p=num2cell(p);
                if nodenumbers
                    s='';
                end
            case 'nodenumbers'
                s='';
                nodenumbers=1;
        end
    end
end

n=length(p);
rr=1;
drr=0.6*(1./(1+exp(2^1.2-n.^1.2))+1)*pi*rr/n;
cr=0.15;
arcr=selfexcitation*0.2*drr;
bc=[1 1 1];

figure
axis([-1 1 -1 1]*(rr+(0.5+cr)*drr+arcr))
axis square
axis off
set(gcf,'Color',bc)

drawnet(p,[0;0],rr,0,s,selfexcitation,grayscale)

function drawnet(p,c,rr,zp,s,selfexcitation,grayscale)
n=length(p);

drr=0.6*(1./(1+exp(2^1.2-n.^1.2))+1)*pi*rr/n;
t=(0:n-1)*2*pi/n;
xy=bsxfun(@plus,rr*[cos(t);sin(t)],c);
bc=[1 1 1];

% generate color vector
sat=[1 0.75 0.5 0.25];
brt=ones(1,n)-grayscale;
cc=hsv2rgb([linspace(0,1-1/n,n);sat(mod(0:n-1,length(sat))+1);brt]');

% draw n units, circles for inhibitory connections, and connection lines
rng=linspace(0,2*pi);
unit=0.5*drr*[cos(rng);sin(rng)];
cr=0.15;
crr=(0.5+cr)*drr;
circ=cr*unit;
z=zeros(1,length(rng))+zp;
arcr=selfexcitation*0.2*drr;
arcrng=linspace(-pi,pi-acos(0.5*arcr/drr));
ang=2*asin(0.5*arcr/crr);
for i=1:n
    h=patch(unit(1,:)+xy(1,i),unit(2,:)+xy(2,i),z,cc(i,:));
    set(h,'EdgeColor',cc(i,:),'LineSmoothing','on')
    if iscell(p{i,i}) 
        set(h,'FaceColor',bc,'LineWidth',1.5)
        rrr=0.9*0.5*drr/(1+(0.5+cr+selfexcitation*0.2)*0.9*pi/length(p{i,i}));
        ss=[s num2str(i)];
        drawnet(p{i,i},xy(:,i),rrr,zp+1,[ss ','],selfexcitation,grayscale)
        %{
        if selfexcitation
            h=text(xy(1,i)+crr*cos(t(i)),xy(2,i)+crr*sin(t(i)),zp+1,ss);
            set(h,'FontUnits','normalized','HorizontalAlignment','center')
            ext=get(h,'Extent');
            set(h,'FontSize',0.4*drr*get(h,'FontSize')/max(ext(3:4)))
        end
        %}
    elseif p{i,i}>0
        if selfexcitation
            xyi=xy(:,i)+crr*[cos(t(i)+ang);sin(t(i)+ang)];
            arcx=xy(1,i)+arcr*cos(t(i)+arcrng)+crr*cos(t(i));
            arcy=xy(2,i)+arcr*sin(t(i)+arcrng)+crr*sin(t(i));
            arcz=0*arcx-0.1+zp;
            line(arcx,arcy,arcz,'Color',cc(i,:),'LineSmoothing','on')
            h=patch(circ(1,:)+xyi(1),circ(2,:)+xyi(2),z,cc(i,:));
            set(h,'EdgeColor',cc(i,:),'LineSmoothing','on')
        end
        if iscell(s)
            ss=s{i};
        else
            ss=[s num2str(i)];
        end
        h=text(xy(1,i),xy(2,i),zp+1,ss,'Color',[0 0 0]+grayscale);
        set(h,'FontUnits','normalized','HorizontalAlignment','center')
        ext=get(h,'Extent');
        set(h,'FontSize',0.8*drr*get(h,'FontSize')/max(ext(3:4)))
    end
	for j=i+1:n
        if p{j,i}>0
            theta=(j+i-2)*pi/n;
            if p{i,j}>0
                xyi=xy(:,j)+crr*[sin(theta-0.25*pi/n);-cos(theta-0.25*pi/n)];
                xyj=xy(:,i)-xyi+xy(:,j);
                xy2=xy(:,i)+crr*[-sin(theta+0.25*pi/n);cos(theta+0.25*pi/n)];
                xxyy=xyi-xy2;
                dd=sqrt(xxyy(1)^2+xxyy(2)^2);
                dx=linspace(0,dd,20)';
                dy=0.2*dx.*(dd-dx)/dd;
                phi=atan2(xxyy(2),xxyy(1));
                xx=[xy(1,i);xy2(1)+cos(phi)*dx-sin(phi)*dy];
                yy=[xy(2,i);xy2(2)+sin(phi)*dx+cos(phi)*dy];
                zz=0*xx-0.1+zp;
                xx(:,2)=xy(1,i)-xx+xy(1,j);
                yy(:,2)=xy(2,i)-yy+xy(2,j);
                h=patch(circ(1,:)+xyj(1),circ(2,:)+xyj(2),z,cc(j,:));
                set(h,'EdgeColor',cc(j,:),'LineSmoothing','on')
            else
                xyi=xy(:,j)+crr*[sin(theta);-cos(theta)];
                xx=[xy(1,i);xyi(1);xyi(1)];
                yy=[xy(2,i);xyi(2);xyi(2)];
                zz=(zp-0.1)+[0;0;0];
            end
            line(xx,yy,zz,'Color',bc,'LineWidth',3,'LineSmoothing','on')
           	h=line(xx,yy,zz,'Color',cc(j,:),'LineSmoothing','on');
            set(h(1),'Color',cc(i,:))
            h=patch(circ(1,:)+xyi(1),circ(2,:)+xyi(2),z,cc(i,:));
            set(h,'EdgeColor',cc(i,:),'LineSmoothing','on')
        elseif p{i,j}>0
            theta=(j+i-2)*pi/n;
            xyj=xy(:,i)+crr*[-sin(theta);cos(theta)];
            xx=[xy(1,j);xyj(1);xyj(1)];
            yy=[xy(2,j);xyj(2);xyj(2)];
            zz=(zp-0.1)+[0;0;0];
            line(xx,yy,zz,'Color',bc,'LineWidth',3,'LineSmoothing','on')
            line(xx,yy,zz,'Color',cc(j,:),'LineSmoothing','on')
            h=patch(circ(1,:)+xyj(1),circ(2,:)+xyj(2),z,cc(j,:));
            set(h,'EdgeColor',cc(j,:),'LineSmoothing','on')
        end
	end
end

function [rho rt]=buildrho(p)
m=length(p);
rho=zeros(m);
rt=cell(1,m);
j=1;
for i=1:m
    rho(j,j+1:end)=cell2mat(p(i,i+1:end));
    rho(j+1:end,j)=cell2mat(p(i+1:end,i));
    if iscell(p{i,i})
        [q qt]=buildrho(p{i,i});
        lq=length(q);
        x=ones(1,lq);
        rho=[rho(1:j-1,1:j-1)	rho(1:j-1,j)*x      rho(1:j-1,j+1:end);
             x'*rho(j,1:j-1)	q                   x'*rho(j,j+1:end);
             rho(j+1:end,1:j-1)	rho(j+1:end,j)*x	rho(j+1:end,j+1:end)];
        for k=1:lq
        	rt{j+k-1}=[num2str(i) ',' qt{k}];
        end
        j=j+lq;
    else
        rho(j,j)=p{i,i};
        rt{j}=num2str(i);
        j=j+1;
    end
end