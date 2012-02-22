function pretty_plot( plt )
% Plot and save a figure as a pdf in a standardized way
% pretty_plot( plt ) will modify the plotted figure and get it ready to be
% saved as a pdf file.
% plt is a structure containing the following elements:
%   plt.fig    = int      % The figure number of the plot
%   plt.legend = {cell}   % Cell list of legend entries(can include option)
%   plt.xlabel = 'char'   % Label for the x-axis
%   plt.ylabel = 'char'   % Label for the y-axis
%   plt.latex  = Logical  % Render text as Latex?
%   plt.pdf    = Logical  % Print to pdf or not
%   plt.file   = 'char'   % File name of the pdf (can include path here)
%   plt.FSn    = int      % FontSize labels/legend
%   plt.FSa    = int      % FontSize axis
%   plt.wide   = float    % Width of the x-axis in the plot
%   plt.AR     = float    % Aspect ratio of the plot
%   plt.zbuf   = float    % Buffer on left/bottom of plot between axis/edge
%   plt.fbuf   = float    % Buffer on right/top of plot between axis/edge
%   plt.misc   = str      % misc. matlab code to be executed before pdf


% Default list of values here
pltD.fig    = 1;      % The figure number of the plot
pltD.pdf    = true;  % Print to pdf or not
pltD.latex  = true;
pltD.file   = 'pretty_plot';   % File name of the pdf (can include path here)
pltD.FSn    = 25;      % FontSize labels
pltD.FSa    = 20;      % FontSize axis
pltD.wide   = 600;     % Width of the x-axis in the plot
pltD.AR     = 4/3;     % Aspect ratio of the plot
pltD.zbuf   = 0.15;    % Buffer on left/bottom of plot between axis/edge
pltD.fbuf   = 0.05;    % Buffer on right/top of plot between axis/edge


% Check the existance of options here
if ( ~ exist('plt') )
    disp('Warning: No options set, default value will be used');
    plt = pltD;
end


%% fig %%
if ( ~ isfield(plt,'fig'))
    plt.fig = pltD.fig;
    disp(['Using default figure value: ' , num2str(plt.fig) ]);
else
    disp(['Using figure value: ', num2str(plt.fig)]);
end

%% legend %%
if ( ~ isfield(plt,'legend'))
    plt.legend_on = false;
    disp('No legend will be used');
else
    plt.legend_on = true;
end

%% xlabel %%
if ( ~ isfield(plt,'xlabel'))
    plt.xlabel_on = false;
    disp('No xlabel will be used');
else
    plt.xlabel_on = true;
end

%% ylabel %%
if ( ~ isfield(plt,'ylabel'))
    plt.ylabel_on = false;
    disp('No ylabel will be used');
else
    plt.ylabel_on = true;
end

%% latex %%
if ( ~ isfield(plt,'latex'))
    plt.latex = pltD.latex;
end

%% pdf %%
if ( ~ isfield(plt,'pdf'))
    plt.pdf = pltD.pdf;
end

%% file %%
if ( ~ isfield(plt,'file'))
    plt.file = pltD.file;
end
disp(['PDF file will be saved as: ',plt.file,'.pdf']);

%% FSn %%
if ( ~ isfield(plt,'FSn'))
    plt.FSn = pltD.FSn;
end

%% FSa %%
if ( ~ isfield(plt,'FSa'))
    plt.FSa = pltD.FSa;
end

%% wide %%
if ( ~ isfield(plt,'wide'))
    plt.wide = pltD.wide;
end

%% AR %%
if ( ~ isfield(plt,'AR'))
    plt.AR = pltD.AR;
end

%% zbuf %%
if ( ~ isfield(plt,'zbuf'))
    plt.zbuf = pltD.zbuf;
end

%% fbuf %%
if ( ~ isfield(plt,'fbuf'))
    plt.fbuf = pltD.fbuf;
end

% Select the desired figure
figure(plt.fig);

% Set the legend
if plt.legend_on
    lstr = '';
    for i=1:size(plt.legend,2)
        if ( i ~= size(plt.legend,2) )
            lstr = [lstr,'''',plt.legend{i},''','];
        else
            lstr = [lstr,'''',plt.legend{i},''''];
        end
    end
    eval( ['h = legend(',lstr,');'] );
    if plt.latex
        set(h,'Interpreter','latex','FontSize',plt.FSn);
    else
        set(h,'FontSize',plt.FSn);
    end
end

% Set the xlabel
if plt.xlabel_on
    h = xlabel( plt.xlabel );
    if plt.latex
        set(h,'Interpreter','latex','FontSize',plt.FSn);
    else
        set(h,'FontSize',plt.FSn);
    end
end

% Set the ylabel
if plt.ylabel_on
    h = ylabel( plt.ylabel );
    if plt.latex
        set(h,'Interpreter','latex','FontSize',plt.FSn);
    else
        set(h,'FontSize',plt.FSn);
    end
end

% Set the axis to the right size
set(gca,'FontSize',plt.FSa);

% Set the plot window to the right size
plt.high = plt.wide / plt.AR;
Nx = plt.wide + (plt.zbuf + plt.fbuf) * plt.wide;
Ny = plt.high + (plt.zbuf + plt.fbuf) * plt.wide;

% Screen Width info
aa = get(0,'ScreenSize');
sW = aa(3);
sH = aa(4);

h = figure(plt.fig);
set(h,'Position',[(sW-Nx)/2,(sH-Ny)/2,Nx,Ny])
%set(h,'PaperSize',[8,6.4]);
set(h,'PaperPositionMode','auto');

% Now set the viewable plot box to be exatly desired size
set(gca,'Position',[plt.zbuf,plt.zbuf*plt.AR,1-(plt.zbuf + plt.fbuf), ...
    1-(plt.zbuf + plt.fbuf)*plt.AR]);

% Misc code execute
if isfield(plt,'misc')
    eval(plt.misc);
end

% Save the plot as a pdf
if plt.pdf
    fname = [ plt.file , '.eps' ];
    print('-depsc2',fname)
    eps2pdf(fname)
    delete(fname)
end


end