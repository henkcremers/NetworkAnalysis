%% color coding for modules and some color maps (cmaps)
% =========================================================================

% https://www.rapidtables.com/web/color/RGB_Color.html
colors = {....,
    [0 100 0],[50 205 50],[144 238 144];...        % 1. green 
    [139 0 0],[255 0 0],[255 99 71];...            % 2. red
    [0 0 139],[0 0 255],[173 216 230];...          % 3. blue
    [255 140 0],[255 170 0],[245 222 179];...      % 4. orange
    [148 0 211],[238 130 238],[218 112 214];...    % 5. purple 
    [105 105 105],[169 169 169],[211 211 211];...  % 6. grey
    [255 20 147],[255 105 180],[255 182 193];...   % 7. pink
    [255 255 0],[255 255 51],[255 255 102];...     % 8. yellow
    };

% rescale codes 
[nr nc] = size(colors);
for j = 1:nr; for i = 1:nc; colors{j,i} = colors{j,i}./256; end; end

%%  plot the colors
% colordat = ones(nc,nr);
% imh = imagesc(colordat);
% figure
% for j = 1:nr;
%     for i = 1:nc;
%         s = scatter(i,nr+1-j,10000,'filled','MarkerFaceColor',colors{j,i})
%         hold on
%     end
% end

%% red-blue color codig 
grade = 10;
n = 64;
mid = (n*0.5);
cmap_hotcold = ones(n,3);
% red scale 
for i = 1:grade;
    
    % first 
    val = i/grade;
    cmap_hotcold(i,1) = 1-val+0.1;
    
    % second 
    cmap_hotcold(i+grade,1) = 0;
    cmap_hotcold(i+grade,2) = 1-val+0.1;
    
    % third 
    cmap_hotcold(i+2*grade,1) = 0;
    cmap_hotcold(i+2*grade,2) = 0;
    cmap_hotcold(i+2*grade,3) = 1-val+0.1;
    
    % black middle
    cmap_hotcold(mid-1,:) = [0.1 0.1 0.1];
    cmap_hotcold(mid,:)   = [0.1 0.1 0.1];
    cmap_hotcold(mid+1,:) = [0.1 0.1 0.1];
    
    % fourth
    cmap_hotcold(33+i,1) = val;
    cmap_hotcold(33+i,2) = 0;
    cmap_hotcold(33+i,3) = 0;
    
    % fifth
    cmap_hotcold(33+i+grade,1) = 1;
    cmap_hotcold(33+i+grade,2) = val;
    cmap_hotcold(33+i+grade,3) = 0;
    
    % sixth
    cmap_hotcold(33+i+2*grade,1) = 1;
    cmap_hotcold(33+i+2*grade,2) = 1;
    cmap_hotcold(33+i+2*grade,3) = val;

end
cmap_hotcold = cmap_hotcold(2:62,:);
% cmpa_hot = cmap_hotcold(,:);

%% red-blue
%if nargin < 1, m = size(get(gcf,'colormap'),1); end
m = 64;
% if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
% else
%     % From [0 0 1] to [1 1 1] to [1 0 0];
%     m1 = floor(m*0.5);
%     r = (0:m1-1)'/max(m1,1);
%     g = r;
%     r = [r; ones(m1+1,1)];
%     g = [g; 1; flipud(g)];
%     b = flipud(r);
% end
cmap_redblue = [r g b]; 
cmap_redblueblack = cmap_redblue;
cmap_red = cmap_redblue(32:end,:);
   
%% markers 
markers = {'o' 's' 'd'};