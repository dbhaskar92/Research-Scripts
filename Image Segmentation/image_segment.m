%{ 
%   Authors: Dhananjay Bhaskar <dbhaskar92@math.ubc.ca>
%   Last modified: Mar, 2018
%   Description: Segment cancer cell line images (Phase contrast microscopy)
%   Tested on MATLAB R2011a
%}

%% Main
function [] = image_segment()

    nfigs = 1;

    test_32(nfigs, 0);
    test_33(nfigs, 0);
    test_34(nfigs, 0);
    test_35(nfigs, 0);
    test_36(nfigs, 0);
    test_37(nfigs, 0);
    test_38(nfigs, 0);
    test_39(nfigs, 0);
    test_40(nfigs, 0);
    test_41(nfigs, 0);
    test_42(nfigs, 0);
    test_43(nfigs, 0);
    test_44(nfigs, 0);
    test_45(nfigs, 0);
    test_46(nfigs, 0);
    % 47 missing
    test_48(nfigs, 0);
    test_49(nfigs, 0);
    test_50(nfigs, 0);
    test_51(nfigs, 0);
    test_52(nfigs, 0);
    % 53 missing
    test_54(nfigs, 0);
    test_55(nfigs, 0);
    test_56(nfigs, 0);
    % 57 missing
    test_58(nfigs, 0);
    test_59(nfigs, 0);
    test_60(nfigs, 0);
    test_61(nfigs, 0);
    test_62(nfigs, 0);
    % 63 missing
    test_64(nfigs, 0);
    % 65 missing
    test_66(nfigs, 0);
    test_67(nfigs, 0);
    % 68, 69 missing
    test_70(nfigs, 0);
    test_71(nfigs, 0);
    test_72(nfigs, 0);
    test_73(nfigs, 0);
    test_74(nfigs, 0);
    test_75(nfigs, 0);
    test_76(nfigs, 0);
    test_77(nfigs, 0);
    test_78(nfigs, 0);
    
    function [nfigs] = test_32(nfigs, display)
        MIAPaCa_32 = strcat('dataset', filesep, 'MIAPaCa_32.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_32, 1.0, 6, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 5
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_32, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_32', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_32_segmented.mat', [6 9 12 14 15 18 19 20 21 23 24 25 27 34], w_cells);
        save_PIF('MIAPaCa_32_segmented.mat', 'MIAPaCa_32.pif');
    end

    function [nfigs] = test_33(nfigs, display)
        MIAPaCa_33 = strcat('dataset', filesep, 'MIAPaCa_33.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_33, 1.0, 7, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 5
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_33, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_33', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_33_segmented.mat', [3 7 10 14 16 18 25], w_cells);
        save_PIF('MIAPaCa_33_segmented.mat', 'MIAPaCa_33.pif');
    end

    function [nfigs] = test_34(nfigs, display)
        MIAPaCa_34 = strcat('dataset', filesep, 'MIAPaCa_34.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_34, 1.0, 12, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 12
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_34, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_34', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_34_segmented.mat', [2 6 12 16 17], w_cells);
        save_PIF('MIAPaCa_34_segmented.mat', 'MIAPaCa_34.pif');
    end

    function [nfigs] = test_35(nfigs, display)
        MIAPaCa_35 = strcat('dataset', filesep, 'MIAPaCa_35.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_35, 0.6, 8, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 10
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_35, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_35', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_35_segmented.mat', [10 13 14 18 19 21], w_cells);
        save_PIF('MIAPaCa_35_segmented.mat', 'MIAPaCa_35.pif');
    end

    function [nfigs] = test_36(nfigs, display)
        MIAPaCa_36 = strcat('dataset', filesep, 'MIAPaCa_36.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_36, 1.0, 10, 4000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 8
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_36, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_36', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_36_segmented.mat', [3 6 8 10 11 13 17 20 22 27 29 30 31 34 36], w_cells);
        save_PIF('MIAPaCa_36_segmented.mat', 'MIAPaCa_36.pif');
    end

    function [nfigs] = test_37(nfigs, display)
        MIAPaCa_37 = strcat('dataset', filesep, 'MIAPaCa_37.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_37, 1.0, 12, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 8
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_37, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_37', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_37_segmented.mat', [2 9 16], w_cells);
        save_PIF('MIAPaCa_37_segmented.mat', 'MIAPaCa_37.pif');
    end

    function [nfigs] = test_38(nfigs, display)
        MIAPaCa_38 = strcat('dataset', filesep, 'MIAPaCa_38.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_38, 1.0, 9, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 6
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_38, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_38', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_38_segmented.mat', [8 11 13 25 26 30], w_cells);
        save_PIF('MIAPaCa_38_segmented.mat', 'MIAPaCa_38.pif');
    end

    function [nfigs] = test_39(nfigs, display)
        MIAPaCa_39 = strcat('dataset', filesep, 'MIAPaCa_39.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_39, 1.0, 10, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 10
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_39, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_39', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_39_segmented.mat', [3 12 13 14 15 20], w_cells);
        save_PIF('MIAPaCa_39_segmented.mat', 'MIAPaCa_39.pif');
    end

    function [nfigs] = test_40(nfigs, display)
        MIAPaCa_40 = strcat('dataset', filesep, 'MIAPaCa_40.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_40, 1.0, 9, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 7
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_40, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_40', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_40_segmented.mat', [2 3 11 12 13 15 18 19 20 26 27 28 29 30 31 34 41 42], w_cells);
        save_PIF('MIAPaCa_40_segmented.mat', 'MIAPaCa_40.pif');
    end

    function [nfigs] = test_41(nfigs, display)
        MIAPaCa_41 = strcat('dataset', filesep, 'MIAPaCa_41.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_41, 1.0, 12, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 12
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_41, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_41', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_41_segmented.mat', [4 5 8 9 11 13 14 18 21 24], w_cells);
        save_PIF('MIAPaCa_41_segmented.mat', 'MIAPaCa_41.pif');
    end

    function [nfigs] = test_42(nfigs, display)
        MIAPaCa_42 = strcat('dataset', filesep, 'MIAPaCa_42.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_42, 1.0, 10, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 8
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_42, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_42', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_42_segmented.mat', [5 7 8 16 18 19 20], w_cells);
        save_PIF('MIAPaCa_42_segmented.mat', 'MIAPaCa_42.pif');
    end

    function [nfigs] = test_43(nfigs, display)
        MIAPaCa_43 = strcat('dataset', filesep, 'MIAPaCa_43.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_43, 1.0, 12, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 9
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_43, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_43', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_43_segmented.mat', [2 3 9 18 21 23 31], w_cells);
        save_PIF('MIAPaCa_43_segmented.mat', 'MIAPaCa_43.pif');
    end

    function [nfigs] = test_44(nfigs, display)
        MIAPaCa_44 = strcat('dataset', filesep, 'MIAPaCa_44.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_44, 1.0, 12, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 9
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_44, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_44', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_44_segmented.mat', [9 10 16 26 29], w_cells);
        save_PIF('MIAPaCa_44_segmented.mat', 'MIAPaCa_44.pif');
    end

    function [nfigs] = test_45(nfigs, display)
        MIAPaCa_45 = strcat('dataset', filesep, 'MIAPaCa_45.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_45, 1.0, 10, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 6
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_45, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_45', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_45_segmented.mat', [15 20 21 23 25], w_cells);
        save_PIF('MIAPaCa_45_segmented.mat', 'MIAPaCa_45.pif');
    end

    function [nfigs] = test_46(nfigs, display)
        MIAPaCa_46 = strcat('dataset', filesep, 'MIAPaCa_46.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_46, 1.0, 14, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 10
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_46, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_46', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_46_segmented.mat', [3 5 11 12 13], w_cells);
        save_PIF('MIAPaCa_46_segmented.mat', 'MIAPaCa_46.pif');
    end

    function [nfigs] = test_48(nfigs, display)
        MIAPaCa_48 = strcat('dataset', filesep, 'MIAPaCa_48.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_48, 1.0, 8, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 6
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_48, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_48', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_48_segmented.mat', [5 9 11 14 28], w_cells);
        save_PIF('MIAPaCa_48_segmented.mat', 'MIAPaCa_48.pif');     
    end

    function [nfigs] = test_49(nfigs, display)
        MIAPaCa_49 = strcat('dataset', filesep, 'MIAPaCa_49.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_49, 1.0, 10, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 10
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_49, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_49', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_49_segmented.mat', [14 16 21 27 29 30 33 34 35], w_cells);
        save_PIF('MIAPaCa_49_segmented.mat', 'MIAPaCa_49.pif');
    end

    function [nfigs] = test_50(nfigs, display)
        MIAPaCa_50 = strcat('dataset', filesep, 'MIAPaCa_50.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_50, 1.0, 10, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 8
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_50, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_50', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_50_segmented.mat', [5 6 8 13 19 21 22 24 25 26], w_cells);
        save_PIF('MIAPaCa_50_segmented.mat', 'MIAPaCa_50.pif');
    end

    function [nfigs] = test_51(nfigs, display)
        MIAPaCa_51 = strcat('dataset', filesep, 'MIAPaCa_51.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_51, 0.7, 12, 3000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 15
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_51, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_51', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_51_segmented.mat', [3 5 6], w_cells);
        save_PIF('MIAPaCa_51_segmented.mat', 'MIAPaCa_51.pif');     
    end

    function [nfigs] = test_52(nfigs, display)
        MIAPaCa_52 = strcat('dataset', filesep, 'MIAPaCa_52.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_52, 1.0, 10, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 6
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_52, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_52', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_52_segmented.mat', [5 6 8 10 13 14 16 18 19 20 23], w_cells);
        save_PIF('MIAPaCa_52_segmented.mat', 'MIAPaCa_52.pif');
    end

    function [nfigs] = test_54(nfigs, display)
        MIAPaCa_54 = strcat('dataset', filesep, 'MIAPaCa_54.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_54, 0.6, 8, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 6
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_54, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_54', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_54_segmented.mat', [14 18 19 22 25 26 29], w_cells);
        save_PIF('MIAPaCa_54_segmented.mat', 'MIAPaCa_54.pif');
    end

    function [nfigs] = test_55(nfigs, display)
        MIAPaCa_55 = strcat('dataset', filesep, 'MIAPaCa_55.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_55, 0.8, 7, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 4
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_55, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_55', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_55_segmented.mat', [9 10 11 12 14 15 17 19 24 25 27], w_cells);
        save_PIF('MIAPaCa_55_segmented.mat', 'MIAPaCa_55.pif'); 
    end

    function [nfigs] = test_56(nfigs, display)
        MIAPaCa_56 = strcat('dataset', filesep, 'MIAPaCa_56.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_56, 0.45, 3, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 2);
        for i = 1 : 10
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_56, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_56', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_56_segmented.mat', [2 9 21 22 26 28 34 40], w_cells);
        save_PIF('MIAPaCa_56_segmented.mat', 'MIAPaCa_56.pif'); 
    end

    function [nfigs] = test_58(nfigs, display)
        MIAPaCa_58 = strcat('dataset', filesep, 'MIAPaCa_58.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_58, 0.45, 3, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 1);
        for i = 1 : 11
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_58, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_58', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_58_segmented.mat', [4 9 20 26 32 35], w_cells);
        save_PIF('MIAPaCa_58_segmented.mat', 'MIAPaCa_58.pif'); 
    end

    function [nfigs] = test_59(nfigs, display)
        MIAPaCa_59 = strcat('dataset', filesep, 'MIAPaCa_59.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_59, 0.7, 7, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 5
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_59, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_59', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_59_segmented.mat', [2 3 10 11 13 19 24 25], w_cells);
        save_PIF('MIAPaCa_59_segmented.mat', 'MIAPaCa_59.pif'); 
    end

    function [nfigs] = test_60(nfigs, display)
        MIAPaCa_60 = strcat('dataset', filesep, 'MIAPaCa_60.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_60, 0.6, 8, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 9
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_60, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_60', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_60_segmented.mat', [14 18 21 25 26 31 37 39 41], w_cells);
        save_PIF('MIAPaCa_60_segmented.mat', 'MIAPaCa_60.pif');
    end

    function [nfigs] = test_61(nfigs, display)
        MIAPaCa_61 = strcat('dataset', filesep, 'MIAPaCa_61.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_61, 0.6, 7, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 2);
        for i = 1 : 11
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_61, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_61', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_61_segmented.mat', [3 4 5 10 12 15 19], w_cells);
        save_PIF('MIAPaCa_61_segmented.mat', 'MIAPaCa_61.pif');
    end

    function [nfigs] = test_62(nfigs, display)
        MIAPaCa_62 = strcat('dataset', filesep, 'MIAPaCa_62.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_62, 1.0, 9, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 2);
        for i = 1 : 9
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_62, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_62', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_62_segmented.mat', [6 11 14 19 20 21 24 26 27], w_cells);
        save_PIF('MIAPaCa_62_segmented.mat', 'MIAPaCa_62.pif');
    end

    function [nfigs] = test_64(nfigs, display)
        MIAPaCa_64 = strcat('dataset', filesep, 'MIAPaCa_64.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_64, 0.8, 6, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 8
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_64, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_64', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_64_segmented.mat', [3 6 8 11 12 13 14], w_cells);
        save_PIF('MIAPaCa_64_segmented.mat', 'MIAPaCa_64.pif');
    end

    function [nfigs] = test_66(nfigs, display)
        MIAPaCa_66 = strcat('dataset', filesep, 'MIAPaCa_66.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_66, 0.8, 8, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 2);
        for i = 1 : 13
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_66, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_66', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_66_segmented.mat', [3 4 6 10 13 19], w_cells);
        save_PIF('MIAPaCa_66_segmented.mat', 'MIAPaCa_66.pif');
    end

    function [nfigs] = test_67(nfigs, display)
        MIAPaCa_67 = strcat('dataset', filesep, 'MIAPaCa_67.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_67, 0.75, 8, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 9
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_67, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_67', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_67_segmented.mat', [7 12 13 21 22 24 25], w_cells);
        save_PIF('MIAPaCa_67_segmented.mat', 'MIAPaCa_67.pif');
    end

    function [nfigs] = test_70(nfigs, display)
        MIAPaCa_70 = strcat('dataset', filesep, 'MIAPaCa_70.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_70, 0.8, 9, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 8
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_70, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_70', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_70_segmented.mat', [2 3 4 6 7 11 18 23 24 26 28], w_cells);
        save_PIF('MIAPaCa_70_segmented.mat', 'MIAPaCa_70.pif');
    end

    function [nfigs] = test_71(nfigs, display)
        MIAPaCa_71 = strcat('dataset', filesep, 'MIAPaCa_71.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_71, 0.8, 8, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 7
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_71, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_71', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_71_segmented.mat', [6 7 11 12 17 22 24], w_cells);
        save_PIF('MIAPaCa_71_segmented.mat', 'MIAPaCa_71.pif');
    end

    function [nfigs] = test_72(nfigs, display)
        MIAPaCa_72 = strcat('dataset', filesep, 'MIAPaCa_72.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_72, 0.8, 9, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 3);
        for i = 1 : 12
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_72, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_72', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_72_segmented.mat', [2 3 8 12 14 16 17], w_cells);
        save_PIF('MIAPaCa_72_segmented.mat', 'MIAPaCa_72.pif');
    end

    function [nfigs] = test_73(nfigs, display)
        MIAPaCa_73 = strcat('dataset', filesep, 'MIAPaCa_73.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_73, 0.8, 8, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 7
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_73, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_73', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_73_segmented.mat', [7 8 9 11 12 13 16 18 21], w_cells);
        save_PIF('MIAPaCa_73_segmented.mat', 'MIAPaCa_73.pif');
    end

    function [nfigs] = test_74(nfigs, display)
        MIAPaCa_74 = strcat('dataset', filesep, 'MIAPaCa_74.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_74, 0.8, 5, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 8
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_74, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_74', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_74_segmented.mat', [4 5 11 15], w_cells);
        save_PIF('MIAPaCa_74_segmented.mat', 'MIAPaCa_74.pif');
    end

    function [nfigs] = test_75(nfigs, display)
        MIAPaCa_75 = strcat('dataset', filesep, 'MIAPaCa_75.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_75, 0.8, 5, 2000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 4
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_75, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_75', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_75_segmented.mat', [3 4 5 7 9 11 12 13 20 21 24 25], w_cells);
        save_PIF('MIAPaCa_75_segmented.mat', 'MIAPaCa_75.pif');
    end

    function [nfigs] = test_76(nfigs, display)
        MIAPaCa_76 = strcat('dataset', filesep, 'MIAPaCa_76.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_76, 0.9, 9, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 7
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_76, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_76', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_76_segmented.mat', [5 10 13 18 19 25], w_cells);
        save_PIF('MIAPaCa_76_segmented.mat', 'MIAPaCa_76.pif');
    end

    function [nfigs] = test_77(nfigs, display)
        MIAPaCa_77 = strcat('dataset', filesep, 'MIAPaCa_77.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_77, 1, 8, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 6
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_77, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_77', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_77_segmented.mat', [3 4 5 10 13 17], w_cells);
        save_PIF('MIAPaCa_77_segmented.mat', 'MIAPaCa_77.pif');
    end

    function [nfigs] = test_78(nfigs, display)
        MIAPaCa_78 = strcat('dataset', filesep, 'MIAPaCa_78.tif');
        % edge detection using mathematical morphology
        [mm_fg, mm_outline, nfigs] = morphological_segment(MIAPaCa_78, 0.8, 9, 1000, nfigs, display);
        % compute foreground markers
        foreground_markers = zeros(size(mm_fg,1), size(mm_fg,2));
        foreground_markers(mm_fg > 0) = 1;
        sedisk = strel('disk', 4);
        for i = 1 : 11
            foreground_markers = imerode(foreground_markers, sedisk);
        end
        foreground_markers = bwareaopen(foreground_markers, 10);
        % marker-based watershed
        [w_cells, w_borders, nfigs] = watershed_segment(MIAPaCa_78, foreground_markers, nfigs, display);
        % plot results
        [nfigs] = plot_results(mm_fg, mm_outline, w_cells, w_borders, 'MIAPaCa_78', nfigs);
        % save correct segmentations
        save_segmentation('MIAPaCa_78_segmented.mat', [15 16 19 21 23 25], w_cells);
        save_PIF('MIAPaCa_78_segmented.mat', 'MIAPaCa_78.pif');
    end
	
end

%% Plot results
function [fig_cnt] = plot_results(mm_fg, mm_outline, w_cells, w_borders, fname, fig_cnt)

    mm_fg_color = label2rgb(mm_fg, 'jet', [.7 .7 .7], 'shuffle');
    mm_outline_color = label2rgb(imdilate(mm_outline, ones(3,3)), 'jet', [.7 .7 .7], 'shuffle');
    w_cells_color = label2rgb(w_cells, 'jet', [.7 .7 .7], 'shuffle');
    w_borders_color = label2rgb(imdilate(w_borders, ones(3,3)), 'jet', [.7 .7 .7], 'shuffle');

    figure(fig_cnt)

    subplot(2,2,1), imshow(mm_fg_color), title('MM Foreground')
    subplot(2,2,2), imshow(w_cells_color), title('Watershed (Marker) Foreground')

    mm_plot = subplot(2,2,3); imshow(mm_outline_color) %, title('MM Outline')
    classes = numel(unique(reshape(mm_fg, [1, numel(mm_fg)])));
    for i = 1 : classes-1
        tmp = zeros(size(mm_fg,1), size(mm_fg,2));
        tmp(mm_fg == i) = 1;
        struct_array = regionprops(tmp, 'centroid');
        if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
            if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
                centroid = cat(1, struct_array.Centroid);
                txt = text(centroid(1), centroid(2), int2str(i));
                set(txt, 'fontsize', 12);
            end
        end
    end

    wm_plot = subplot(2,2,4); imshow(w_borders_color)
    classes = numel(unique(reshape(w_cells, [1, numel(w_cells)])));

    for i = 1 : classes
        tmp = zeros(size(w_cells,1), size(w_cells,2));
        tmp(w_cells == i) = 1;
        struct_array = regionprops(tmp, 'centroid');
        if ~isempty(struct_array) && isfield(struct_array, 'Centroid')
            if ~isempty(struct_array.Centroid) && numel(extractfield(struct_array, 'Centroid')) == 2
                centroid = cat(1, struct_array.Centroid);
                txt = text(centroid(1), centroid(2), int2str(i));
                set(txt, 'fontsize', 12);
            end
        end
    end

    % Print outlines

    fig = figure;
    set(mm_plot, 'Visible', 'off');
    mm_fig = copyobj(mm_plot, fig);
    set(mm_fig, 'Position', get(0, 'DefaultAxesPosition'));
    saveas(mm_fig, strcat(fname, '_MM', '.png'));

    fig = figure;
    set(wm_plot, 'Visible', 'off');
    wm_fig = copyobj(wm_plot, fig);
    set(wm_fig, 'Position', get(0, 'DefaultAxesPosition'));
    saveas(wm_fig, strcat(fname, '_WS', '.png'));

    fig_cnt = fig_cnt + 1;
    
end

%% Save segmentation
function save_segmentation(filename, cell_labels, labelled_image)

    segmentations = zeros(size(labelled_image,1), size(labelled_image,2));
    max_cid = 0;

    if exist(filename, 'file')
        load(filename, 'segmentations', 'max_cid');
        assert(isequal(size(segmentations),size(labelled_image)), 'Matrix dimensions do not match.');  
    end

    num_correct_segmentations = numel(cell_labels);

    for i = 1 : num_correct_segmentations
        max_cid = max_cid + 1;
        class = cell_labels(i);
        segmentations(labelled_image == class) = max_cid;
    end

    save(filename, 'segmentations', 'max_cid');
    
end

%%  Convert MAT to PIF file
function save_PIF(filename, PIFname)

    segmentations = [];

    if exist(filename, 'file')
        load(filename, 'segmentations');
    else
        print('Error: MAT file does not exist.'); 
    end

    cell_ids = unique(nonzeros(segmentations));
    numPIFrows = nnz(segmentations);

    PIF_data = zeros(numPIFrows, 5);
    ind = 1;
    for i = 1 : length(cell_ids)
        [rows, cols] = find(segmentations == cell_ids(i));
        for cnt = 1 : length(rows)
            PIF_data(ind, 1) = cell_ids(i);
            PIF_data(ind, 2) = rows(cnt);
            PIF_data(ind, 3) = rows(cnt);
            PIF_data(ind, 4) = cols(cnt);
            PIF_data(ind, 5) = cols(cnt);
            ind = ind + 1;
        end
    end
    
    if exist(PIFname, 'file')
        print('Error: PIF file already exists.');
    else
        fileID = fopen(PIFname, 'w');
        cnt = 1;
        while cnt < ind
            fprintf(fileID, '%d CellU %d %d %d %d\n', PIF_data(cnt,:));
            cnt = cnt + 1;
        end
        fclose(fileID);
    end

end