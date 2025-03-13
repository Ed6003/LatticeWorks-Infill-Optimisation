cd("C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks\TopOpt") %(*@\draftcomment{the prefix for every sub-heading in this script is "opt"}@*)
savePath = 'D:\TechnicalReport\Variable-TPMS-Figures\TopOpt';
% This work is based on DEMO 13/14

load('sample.mat');

model = createpde("structural","static-solid");

% Import geometry from a STEP file
importGeometry(model, 'Planet-18-Teeth-1.8mm-Mod.step');
scale(model.Geometry,1e3); % to mm
material_yield_strength = 50e6; % MPa
material_E = 1666e6; % Pa
material_v = 0.38;

% Visualise the imported geometry
cFigure;
pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.5);
title('Imported Gear Geometry');
xlim([-20 20])
ylim([-20 20])
zlim([0 20])

teeth_faces = [14,24,31,38,48,55,62,72,79,86,96,103,110,120,127,134,144,151];

teeth_distance = 15e-3; % mm
torque = 10; % J (Nm)
teeth_force = (torque / teeth_distance) / numel(teeth_faces); % N

area_mm2 = 66.441 * 1e-6;
pressure = teeth_force / area_mm2;

% Apply pressure load to each gear tooth face
structuralBoundaryLoad(model,"Face",teeth_faces,"Pressure",pressure);

structuralProperties(model,"YoungsModulus",material_E,"PoissonsRatio",material_v);
structuralBC(model,"Face",8,"Constraint","fixed");
generateMesh(model, 'GeometricOrder', 'quadratic','Hmax',0.0005);
fprintf('%g k Elements created',size(model.Mesh.Elements,2)./1000);

figure
pdeplot3D(model)
title("Mesh with Quadratic Tetrahedral Elements");

result = solve(model);

maxVonMises = max(result.VonMisesStress); % As the load is distributed using surface traction and surface traction = F / A
fprintf("Max Von Mises is %g MPa.", (maxVonMises / 1e6) );
fprintf("Safety Factor is %.2f.", (material_yield_strength / (maxVonMises)));

% Data to visualise
meshData = result.Mesh;
nodalData = result.VonMisesStress;
deformationData = result.Displacement;

% Create PDE result visualisation
resultViz3 = pdeviz(meshData,nodalData, ...
    "DeformationData",deformationData, ...
    "DeformationScaleFactor",0, ...
    "ColorLimits",[1700 1386000]);

% Clear temporary variables
clearvars meshData nodalData deformationData

%(*@\codesubsection{3D Grid Interpolation}{opt-3D-grid-interpolation}@*)
x = result.Mesh.Nodes(1, :);
y = result.Mesh.Nodes(2, :);
z = result.Mesh.Nodes(3, :);
vonMisesStress = result.VonMisesStress;

n = 200;
xq = linspace(min(x), max(x), n);
yq = linspace(min(y), max(y), n);
dx = (max(x)-min(x))/(n-1);

% Estimate z for square voxels
nz = round((max(z) - min(z))/dx) + 1;
zq = linspace(min(z), max(z), nz);

[Xq, Yq, Zq] = meshgrid(xq, yq, zq);

% Create an interpolant function (F)
F = scatteredInterpolant(x', y', z', vonMisesStress, 'linear', 'none');

% Evaluate interpolant at the meshgrid
stressGrid = F(Xq, Yq, Zq);

stressGrid(isnan(stressGrid)) = 0;

%(*@\codesubsection{Masking}{opt-masking}@*)
alphaVal = 1; % Pa threshold (arbitrary)
shp = alphaShape(x', y', z', alphaVal);
insideMask = inShape(shp, Xq, Yq, Zq);
stressGrid(~insideMask) = NaN;

%(*@\codesubsection{Interpolated Von Mises}{opt-interpolated-von-mises}@*)
cFigure;
surf(Xq(:,:,1), Yq(:,:,1), stressGrid(:,:,1), 'EdgeColor', 'none');
colorbar;
% title('Interpolated Von Mises Stress Distribution at Z = min(Z)');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Von Mises Stress');
view(2);

% slice planes
xslice = [min(Xq(:)), mean(Xq(:)), max(Xq(:))];
yslice = [min(Yq(:)), mean(Yq(:)), max(Yq(:))];
zslice = [min(Zq(:)), mean(Zq(:)), max(Zq(:))];

% slice plot
cFigure;
slice(Xq, Yq, Zq, stressGrid, xslice, yslice, zslice);
shading interp;
colorbar;
% title('3D Slices of Von Mises Stress Distribution');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');

%(*@\codesubsection{Stress Based Infill Distribution}{opt-stress-based-infill-distribution}@*)
p95Stress = prctile(stressGrid(:), 95); % 95th percentile stress
clampedStress = min(stressGrid, p95Stress);
minStress = min(clampedStress(:), [], 'omitnan'); % omit nans just in case

% Normalize using clamped min val
rhoGrid = (clampedStress - minStress) ./ (p95Stress - minStress);

% Clamping check
if any(rhoGrid > 1, 'all') || any(rhoGrid < 0, 'all')
    error('Infill outside 0-100%% range created.');
end

rhoGrid(~insideMask) = NaN; % previously was rhoGrid(isnan(stressGrid)) = NaN;

%(*@\codesubsection{Infill Distribution Plot}{opt-infill-distribution-plot}@*)
% Infill Distribution
cFigure;
surf(Xq(:,:,1), Yq(:,:,1), rhoGrid(:,:,1), 'EdgeColor', 'none');
colorbar;
% title('Infill Distribution at Z = min(Z)');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Von Mises Stress');
view(2);

% Infill Distribution
halfwayZ = ceil(size(Xq,3)/2);
cFigure;
surf(Xq(:,:,halfwayZ), Yq(:,:,halfwayZ), rhoGrid(:,:,halfwayZ), 'EdgeColor', 'none');
colorbar;
% title('Infill Distribution at centroid');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Von Mises Stress');
view(2);

%(*@\codesubsection{Create Solid Gear Teeth}{opt-create-solid-gear-teeth}@*)
radius = 13.95; % minimum radius to gear teeth

% Logical mask of points outside circle
circleMask = (Xq.^2 + Yq.^2) <= radius^2;
rhoGrid(~circleMask & ~isnan(rhoGrid)) = 1; % Set solid infill teeth

cFigure;
surf(Xq(:,:,1), Yq(:,:,1), rhoGrid(:,:,1), 'EdgeColor', 'none');
colorbar;
% title('Infill Distribution at Z = min(Z)');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Von Mises Stress');
view(2);

%(*@\codesubsection{Patch Model Data}{opt-patch-model-data}@*)
v = 1;

% Patch rho information to Elements and Faces
[E,V]=im2patch(rhoGrid,true(size(rhoGrid)),'h',v);
[F,rho_F] = element2patch(E,rhoGrid(:));

% Boundary faces
indBoundaryFaces = tesBoundary(F);
Fb = F(indBoundaryFaces,:);

% define face boundary rho values
if size(Fb,2)==4
    Fb = [Fb(:,[1 2 3]); Fb(:,[3 4 1]);]; % force to be triangles to be compatible with tet volume based in/out testing
    rho_Fb = [rho_F(indBoundaryFaces); rho_F(indBoundaryFaces)];
else
    rho_Fb = rho_F(indBoundaryFaces);
end

%(*@\codesubsection{Create \figurednt{TPMS}}{topopt-create-TPMS}@*)
k = pi; % k is grid size

% Generate TPMS function S
S = sin(k*Xq).*cos(k*Yq) + sin(k*Yq).*cos(k*Zq) + cos(k*Xq).*sin(k*Zq);
S(~insideMask) = 0;

rho_VG = rhoGrid;

% Mapping the density to the equivalent gyroid levelSet field
S_flat = S(:);
sorted_S = sort(S_flat);
n = numel(sorted_S);
percentiles = (1:n)/n; % Percentiles from 0 to 1 (exclusive)

target_percentiles = 1 - rho_VG;
beta = interp1(percentiles, sorted_S, target_percentiles, 'pchip');
beta = max(min(beta, max(sorted_S)), min(sorted_S));

% Generate Sg and extract isosurface
Sg = S - beta;

[Fg,Vg] = isosurface(Xq,Yq,Zq,Sg,0);
[Fgc,Vgc] = isocaps(Xq,Yq,Zq,Sg,0);

[Fsn,Vsn,Csn]=joinElementSets({Fg,Fgc},{Vg,Vgc});
[Fsn,Vsn]=patchCleanUnused(Fsn,Vsn); % Remove unused nodes
[Fsn,Vsn]=mergeVertices(Fsn,Vsn); % Merge nodes

S(~insideMask) = NaN;
sv3(S); colormap warmcold;
Sg(~insideMask) = NaN;
sv3(Sg); colormap warmcold;

%(*@\codesubsection{Create Generated Gear Geometry}{opt-generated-gear-geometry}@*)
cFigure;
patch('Faces', Fsn, 'Vertices', Vsn, 'FaceColor', 'blue', 'EdgeColor', 'none');
view(3); axis equal; camlight; lighting gouraud;
% title('Variable-Density TPMS Structure');

cFigure;
patch('Faces', Fsn, 'Vertices', Vsn, 'FaceColor', 'blue', 'EdgeColor', 'none');
view(2); axis equal; camlight; lighting gouraud;
% title('Variable-Density TPMS Structure');

TR = triangulation(Fsn, Vsn);
stlwrite(TR, fullfile(savePath,'OutputGear.stl'));

cFigure;
scatter3(Xq(:),Yq(:),Zq(:),20,S(:),'filled');
axis tight; axis equal;
colorbar;

cFigure;
scatter3(Xq(:),Yq(:),Zq(:),20,Sg(:),'filled');
axis tight; axis equal;
colorbar;

%(*@\codesubsection{Binary Mask Based on Sign}{topopt-binary-mask-based-on-sign}@*)
% Create binary mask (1 where Sg > 0, 0 where Sg <= 0)
binary_mask = Sg > 0;

x_slice = 0;
y_slice = 0;
z_slice = max(Zq(:))/2;

slice(Xq, Yq, Zq, double(binary_mask), x_slice, y_slice, z_slice);
colormap([0.3010, 0.7450, 0.9330; 0, 0.4470, 0.7410]); % [Blue for 0 (Sg <= 0), Red for 1 (Sg > 0)]
shading flat;

axis tight equal;
colorbar('off'); % Remove colorbar
% title('Sg > 0 (Blue) vs. Sg < 0 (Cyan)');

%(*@\codesubsection{Sliced Mesh Visualisation}{opt-sliced-mesh-visualisation}@*)
meanZ = mean(Vsn(:,3));

P = [0, 0, meanZ]; % Cut plane point
n = [0, 0, 1]; % Z Normal

snapTolerance = mean(patchEdgeLengths(Fsn, Vsn)) / 100;
[Fc, Vc, ~, logicSide, Eb] = triSurfSlice(Fsn, Vsn, [], P, n, snapTolerance);

cFigure; hold on;
gpatch(Fsn, Vsn, 'kw', 'none', 0.01); % transparent original
gpatch(Fc(logicSide,:), Vc, 'rw', 'k', 1); % sliced mesh
gpatch(Eb, Vc, 'none', 'b', 1); % intersections

% Enhance visualisation
axisGeom;
camlight headlight;
colormap gjet;
icolorbar;
% title('Sliced Mesh Visualisation');
drawnow;

%(*@\codesubsection{Re-Save New Figures}{opt-resave-new-figures}@*)
saveFigures(savePath,true,0);