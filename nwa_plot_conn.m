function [] = nwa_plot_conn(C,varargin)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

newfigure = 0;

nn = size(C);
if ~isequal(nn(1),nn(2)); error('not a square matrix'); end
nn = nn(1);

% defaults
method = 'corrplot';
t = 'Correlation Matrix';
for j = 1:nn;
    varnames{j} = ['variable' num2str(j)];
end


% user input
%-----------
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'varnames', varnames = varargin{i+1};
            case 'title', t = varargin{i+1};
            case 'method', method = varargin{i+1};
            case 'modsort'; 
                [Ci Q] = modularity_und(C); 
                disp(['Modularity is: ' num2str(Q)])
                [dump loc] = sort(Ci);
                C = C(loc,loc);
            case 'figure', newfigure = 1;    
        end
    end
end
if newfigure==1;
    figure
end

%% plot method
switch method
    
    case 'corrplot'

        % figure;
        imagesc(1:nn,1:nn,C);
        %load('nwa_colormap.mat')
        nwa_colors;
        % colormap(cmap_hotcold)
        colormap(cmap_redblue)
        %colormap(mycmap(32:end,:))
        % caval = 1;
        caval = max(max(abs(C)));
        caxis([-caval,caval])
        
        if nn<20
            % plot the matrix values
            for i = 1:nn;
                for j = 1:nn;
                    val = round(C(i,j),2);
                    fcol = 'k';
                    if val<0.2; fcol = 'w'; end
                    text(i, j, num2str(val), 'FontSize', 12, 'Color', fcol);
                end;
            end
            
            % fix the axis;
            set(gca, 'XTick', 1:nn, 'XTickLabel', varnames);
            set(gca, 'YTick', 1:nn, 'YTickLabel', varnames);
            
        else 
            axis off
        end
        
        % add title
        title(t)
        
    case 'network'
        
        % thrshold
        C = nwa_proc_conn(C,'threshr',0.05);
        
        % set identity to zero
        C0 = C;
        C0(logical(eye(nn))) = 0;
        [ind1,ind2]=ind2sub(size(C0),find(C0~=0));
        
        % set up and plot the nodes
        theta=linspace(0,2*pi,nn+1);theta=theta(1:end-1);
        [x,y]=pol2cart(theta,1);
        plot(x,y,'.k','markersize',10);hold on
        
        % label the nodes
        txt = cellstr(num2str((1:nn)','%02d'));
        h = text(x.*1.05, y.*1.05, txt, 'FontSize',8);
        % set(h, {'Rotation'},num2cell(theta*180/pi))
        
        % draw the edges.
        for j = 1:length(ind1)
            eX = [x(ind1(j)) x(ind2(j))];
            eY = [y(ind1(j)) y(ind2(j))];
            eW = C(ind1(j),ind2(j));
            if eW > 0; eC = 'g'; else; eC = 'r'; end
            eW = abs(eW);
            plot(eX,eY,eC,'LineWidth',eW);
        end
        
        axis off
        title(t)
        
        %         arrayfun(@(p,q)line([x(p),x(q)],[y(p),y(q)]),ind1,ind2);
        %         axis equal off
        
    case 'hist'
        
        cdat = nwa_reshape(C,'mat2vec');
        dens = sum(cdat~=0)/length(cdat);
        b = [min(cdat):0.005:max(cdat)];
        [h x] = histc(cdat,b);
        [h2 x2] = histc(cdat(cdat~=0),b);
        %h2 = h2.*(max(h)/max(h2));
        subplot(2,1,1)
        plotyy(b,h,b,h2)
%         p = plot(b,h,'k'); %hold on
% 
%         p = plot(b,h2,'r');
%        legend(['all: density:' num2str(dens)],'non-zero')
        
        title('Edge distribution')
        subplot(2,1,2)
        deg = degrees_und(C);
        range = [min(deg):max(deg)];
        [h x] = histc(deg,range);
        p = plot(range,h,'k');
        title('Degree distribution')
           
end

