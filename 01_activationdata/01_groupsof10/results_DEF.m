%--------------------------------------------------------------------------
% Sample Matlab code for creating images with hue and alpha color-mapping.
%
%  Notes for using this code:
%  You must have OpenGL available on your system to use transparency (alpha). 
%  When rendering transparency MATLAB automatically uses OpenGL if it is
%  available. If it is not available, transparency will not display. 
%  See the figure property RendererMode for more information.
%
% EA Allen August 30, 2011
% eallen@mrn.org
%--------------------------------------------------------------------------
SLICES = [26,27,37];
%for i = SLICES
 i=42   
    
    %% 1. Load the AOD_data.mat file with sample data from the fMRI AOD experiment
    %--------------------------------------------------------------------------
    addpath('./Onderzoek/doctoraat/scripts/01_MATLAB/dualcodeExample')
    addpath('./Onderzoek/doctoraat/reports_abstracts/Paper_finaal/03_manuscript/00_pix/')
    fname=strcat('CON_emo_',num2str(i),'.mat');
    %fname=strcat('CON_emoS6_',num2str(i),'.mat');
    load(fname) 
    % For a single axial slice (Z = 2 mm) of data, you should have:
    % Bmap_N_S: 'Difference between Novel and Standard betas averaged over 28 subjects'
    % Tmap_N_S: 'T-statistics for the paired t-test comparing Novel and Standard betas'
    % Pmap_N_S: 'Binary map indicating significance at P<0.001 (fdr corrected)'
    % Underlay: 'Structural image ch2bet from MRIcron, warped to functional data'
    %--------------------------------------------------------------------------

    %% 2. Set some defaults that will affect the appearance of the image
    %--------------------------------------------------------------------------
    % Set the Min/Max values for hue coding
    absmax = max(max(Z)); 
    absmin = min(min(Z)); 
    H_range = [absmin absmax]; % The colormap is symmetric around zero
    % Set the Min/Max T-values for alpha coding
    A_range = [0 1.0];
               
    % Set the labels for the colorbar
    hue_label = 'Z';
    alpha_label = 'reselection rate';

    % Choose a colormap for the underlay
    CM_under = gray(256);

    % Choose a colormap for the overlay
    CM_over = jet(256);
    
    
    % Define selected cluster with min integer
    Bsel =B;
    Bstab=B;
    Bsel(Bsel<99)=0;
    Bstab(Bstab>99)=0;
    Bstab(Bstab<97)=0;
    
    %--------------------------------------------------------------------------

    %% 3. Do the actual plotting
    %--------------------------------------------------------------------------
    % Make a figure and set of axes
    F = figure('Name', 'slice', 'Color', 'k', 'Units', 'Normalized', 'Position', [0.3, 0.4, 0.2, 0.35] ); 
    axes('Position', [0 0 1 1]); 

    % Transform the underlay and beta map to RGB values, based on specified colormaps
    % See function convert_to_RGB() for more information
    U_RGB = convert_to_RGB(A, CM_under);

    statmap = Z;
    statmap(statmap > H_range(2)) = H_range(2);
    %statmap(statmap < 3.09) = 3.08;

    O_RGB = convert_to_RGB(statmap, CM_over);

    % Plot the underlay
    layer1 = image(U_RGB); axis image
    hold on;
    % Now, add the Beta difference map as an overlay
    layer2 = image(O_RGB); axis image

    % Use the T-statistics to create an alpha map (which must be in [0,1])
    alphamap = abs(S);
    alphamap(alphamap > A_range(2)) = A_range(2);
    alphamap(alphamap < A_range(1)) = 0;
    alphamap = alphamap/A_range(2);

    % Adjust the alpha values of the overlay 
    set(layer2, 'alphaData', alphamap);

    % Add some (black) contours to annotate nominal significance
    hold on;
    [C, CH] = contour(Bsel, 1, 'w','LineWidth',1.5);
   % hold on;
    [C2, CH2] = contour(Bstab, 1, 'b','LineWidth',1.5);

    set(gcf,'PaperPositionMode','auto')
    fname2=strcat('./Onderzoek/doctoraat/reports_abstracts/Paper_finaal/03_manuscript/00_pix/brain_',num2str(i), '_', num2str(j*100));
    print(fname2 ,'-dpng','-r0')
    %--------------------------------------------------------------------------

    %% 4. Create a 2D colorbar for the dual-coded overlay
    %--------------------------------------------------------------------------
    G = figure('color', 'k', 'Units', 'Normalized', 'Position', [0.5, 0.4, 0.06, 0.35]);
    x = linspace(A_range(1), A_range(2), 256); 
    % x represents the range in alpha (abs(t-stats))
    y = linspace(H_range(1), H_range(2), size(CM_over,1));
    % y represents the range in hue (beta weight difference)
    [X,Y] = meshgrid(x,y); % Transform into a 2D matrix
    imagesc(x,y,Y); axis xy; % Plot the colorbar
    set(gca, 'Xcolor', 'w', 'Ycolor', 'w')
    colormap(CM_over); 
    alpha(X);
    alpha('scaled');  
    xlabel(alpha_label)
    set(gca, 'YAxisLocation', 'right')
    ylabel(hue_label)
    set(gcf,'PaperPositionMode','auto')
    fname2=strcat('./Onderzoek/doctoraat/reports_abstracts/Paper_finaal/03_manuscript/00_pix/label',num2str(i), '_', num2str(j*100));
    print(fname2 ,'-dpng','-r0')
    %------------------------------------------------------------------------
%end
