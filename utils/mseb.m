function H = mseb(x,y,errBar,lineProps,transparent)
% H = MSEB(x,y,errBar,lineProps,transparent)
%
% Multiple Shaded Error Bars (MSEB) makes a 2-d plot containing multiple 
% lines with pretty shaded error bars.
% This is an extension of the popular shadedErrorBar, by Rob Campbell,
% enabling plotting of multiple data lines with overlapping errorbars, as 
% well as turning off legends for all elements but the main lines. MSEB is
% directed at using the default renderer instead of openGL, which is known
% to cause problems, e.g., with logarithmic axes and not being able to save
% figures as vector graphics in the eps-format.
% To avoid the error bar patches concealing previously plotted main lines, 
% the different elements are plotted in a suitable order. The default
% setting plot edges of overshadowed patches in an non-obtrusive manner so 
% all error bars can be seen but avoids cluttering of the plot (see the 
% third example on how to customise this). The patches are plotted in
% reverse order to make sure the first entry is "on top".
%
% Andreas Trier Poulsen, atpo@dtu.dk
% Technical University of Denmark, Cognitive systems - January 2015
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of C x N, where C is the number of
%     lines to be plotted and N is the number of samples in each line and
%     should be equal to length(X)
% errBar - if a vector we draw symmetric errorbars. If it has a
%          size of [C,length(x),2] then asymmetric error bars are drawn,
%          with row 1 being the upper bar and row 2 being the lower
%          bar. In the present version errBar does not support two function
%          handles.
% lineProps  - [optional. Can also be set as "[]" for default settings].
%              Struct containing fields that define lineproperties for the 
%              plot function. It is possible to only define some of the
%              fields.
%       .col - cell array, where each element defines the colour of each
%              line. This can be done using either strings or three-element
%              RGB vectors e.g., either 'b' or [0 0 1] for blue.
%     .style - linestyle of the lines from y-data. Default is '-'.
%     .width - linewidth of the lines from y-data. Default is 2.
% .edgestyle - linestyle of edges that are overlapped by errorbars from
%              other lines. 
%             
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image. openGL does not support having logarithmic axes.
%
% Outputs
% H - structure with an element for each line entry containing handles to 
%     the generated plot objects (e.g. H(c) contains the handles to the
%     c'th line entry. 

% Error checking
narginchk(3,5)

% Cheking the y data
[C, N] = size(y);
if N==1
    C = 1;
    N = length(y);
    y = y';
end

% Cheking the x data
if isempty(x)
    x=repmat(1:N,[C,1]);
elseif length(x(:))==N && C>1
    x=repmat(x(:)',[C,1]);
end

if (size(x,1) ~= C) || (size(x,2)~=N)
    error('inputs x and y do not have same dimensions')
end

% Checking ErrBar dimensions. If only one error bar is specified then we
% will mirror it, turning it into both upper and lower bars.
if (size(errBar,1) ~= C) || (size(errBar,2)~=N)
    if size(errBar,1)==N && size(errBar,2)==2 && size(errBar,3)~=1 && C==1
        errBar = permute(errBar,[3,1,2]);
    elseif size(errBar,1)==N && size(errBar,2)==1 && size(errBar,3)~=1 && C==1
        errBar=repmat(errBar,[1,2]);
        errBar = permute(errBar,[3,1,2]);
    else
        error('inputs errBar and y do not have same dimensions')
    end
end
if size(errBar,3)==1
    errBar = repmat(errBar,[1,1,2]);
end

% Setting default line options
if nargin<4
    lineProps = [];
end
if nargin<5 || ~isnumeric(transparent)
    transparent=0;
end
if ~isfield(lineProps,'col')
    for c = 1:C
        colours = 'brkgmcy';
        col_ind = rem(c-1,length(colours))+1;
        lineProps.col{c} = colours(col_ind);
    end
end
if ~isfield(lineProps,'style'),	lineProps.style = '-'; end
if ~isfield(lineProps,'width'),	lineProps.width = 2; end
if ~isfield(lineProps,'edgestyle'),	lineProps.edgestyle = '--'; end

% First loop over the number of lines to be plotted
holdStatus=ishold;
if ~holdStatus, hold on,  end

for c = C:-1:1
    H(c).mainLine=scatter(x(c,:),y(c,:), 10,'filled','MarkerFaceColor',lineProps.col{c},'MarkerEdgeColor','none');

    col=get(H(c).mainLine,'MarkerFaceColor');
    patchSaturation=0.20;
    if transparent
        faceAlpha=patchSaturation;
        patchColor = col;
        set(gcf,'renderer','openGL')
    else
        faceAlpha=1;
        patchColor=col+(1-col)*(1-patchSaturation);
        set(gcf,'renderer','painters')
    end
    edgeAlpha = faceAlpha;
    edgeColor = patchColor;
    
    uE=y(c,:)+errBar(c,:,1);
    lE=y(c,:)-errBar(c,:,2);
    
    yP=[lE,fliplr(uE)];
    xP=[x(c,:),fliplr(x(c,:))];
    
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];
    
    H(c).patch=patch(xP,yP,1,"FaceColor",patchColor,"FaceAlpha",faceAlpha,"EdgeColor",edgeColor,"EdgeAlpha",edgeAlpha);
end

for c = 1:C
    delete(H(c).mainLine)
    H(c).mainLine=scatter(x(c,:),y(c,:), 10,'filled','MarkerFaceColor',lineProps.col{c},'MarkerEdgeColor','none');

    set(get(get(H(c).patch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    set(get(get(H(c).mainLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','on')
end
if ~holdStatus, hold off, end
