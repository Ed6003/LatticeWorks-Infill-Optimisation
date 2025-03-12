defaultFolder = fullfile('D:','TechnicalReport','LatticeProperties_LinearFinal');

%(*@ Previously there appears to be a mistake in the Young's modulus
% calculation, since the absolute of th strain was used to plot the young's
% modulus the difference in the strain was halved, ie. if you have +-0.2,
% the range goes from 0.4 to now being 0.2 so it needs to be multiplied by 2
% it also needs to be multiplied by 10 so in total divided by 5 %@*)

original_points = linspace(1.3, -1.3, 64);

exclusions = {'f', 'v', 'c', 'F', 'V', 'meshOutput', 'abaqusData', 'E_effectiveStrain', 'E_effectiveStress'};

for i = linspace(1.3,-1.3,64)

    savePath=fullfile(defaultFolder,sprintf('%.5g_Lattice_Density',i));
    matlabPath=fullfile(savePath,'simulation_results.mat');

    load(matlabPath);

    resultStruct.youngs_modulus.youngs_modulus_xz_mean = resultStruct.youngs_modulus.youngs_modulus_xz_mean/5;
    resultStruct.youngs_modulus.youngs_modulus_xz_median = resultStruct.youngs_modulus.youngs_modulus_xz_median/5;
    resultStruct.youngs_modulus.youngs_modulus_xy_mean = resultStruct.youngs_modulus.youngs_modulus_xy_mean/5;
    resultStruct.youngs_modulus.youngs_modulus_xy_median = resultStruct.youngs_modulus.youngs_modulus_xy_median/5;

    save(matlabPath, 'resultStruct');
    disp(i)
end