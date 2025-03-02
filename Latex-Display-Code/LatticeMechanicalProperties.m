clear; close all; clc; %(*@\draftcomment{the prefix for every sub-heading in this script is "lattice"}@*)

updateSim = true; % will not recompute abaqus simulation, useful for testing fig generation
overwrite = false; % overwrite already simulated files (forces script completion)
coarseMesh = true; % coarse mesh for testing (quick sim)
figureTitles = true;
loadSimDataFromMat = false; % load data from mat and not dat file if sim exists (much faster)
doSaveFigures = true; % skip figure saving and just close them

if ~updateSim && overwrite
    warning('updateSim is false but overwrite is true, figures will be updated but Abaqus simulation will not')
end

cd('C:\Users\rusco\OneDrive - University of Warwick\Admin\Archives\Documents\GitHub\LatticeWorks');

try
    % Start the Data Collector Set
    system('logman start CPU_Log');

    %(*@\codesubsection{Plot Settings}{lattice-plot-settings}@*)
    cMap=[0.6*ones(256,2), linspace(1, 0, 256)'];
    fontSize=15;
    markerSize=20;

    % Define default colours
    blue   = [0, 0.4470, 0.7410];  % "#0072BD"
    cyan  = [0.3010, 0.7450, 0.9330];  % "#4DBEEE"

    %(*@\codesubsection{Control Parameters}{lattice-control-parameters}@*)
    pointSpacing = 0.15; % changing this to see if meshes within 250k

    % Generating Geometry
    for i = linspace(1.3,-1.3,64) % 64

        nFigures = 0;

        levelset = i; % Isosurface level, corresponding to volume fraction

        x_length = 30; y_length = 10; z_length = 10; % geometry
        inputStruct.L=[x_length y_length z_length]; % characteristic length default: inputStruct.L=[3 1 1]
        inputStruct.Ns=100; % number of sampling points, resolution
        inputStruct.surfaceCase='g'; %Surface type
        inputStruct.numPeriods=[9 3 3]; %Number of periods in each direction
        inputStruct.gradType='levelSet';
        inputStruct.GF=[1, 1]; % Gradient factors (constant gradient)

        %(*@\codesubsection{Abaqus Simulations}{lattice-abaqus-simulations}@*)
        resultStruct = struct();

        defaultFolder = 'D:'; % Changed to external hard drive
        savePath=fullfile(defaultFolder,'TechnicalReport','LatticeProperties',sprintf('%.5g_Lattice_Density',i));

        % exits if already completed
        if simulationCompleted(savePath) && ~overwrite
            fprintf('[%.5g] Lattice density already completed, skipping ...\n', i);
            % error('Enabled overwrite required');
            continue;  % Skip to the next iteration
        end

        % ensures folders exist
        [status, msg, msgID] = mkdir(savePath);
        if ~status
            error('Failed to create directory: %s', msg);
        end

        % defining file names
        abaqusInpFileNamePart='Lattice_FEA';
        abaqusInpFileName=fullfile(savePath,[abaqusInpFileNamePart,'.inp']); %INP file name
        abaqusDATFileName=fullfile(savePath,[abaqusInpFileNamePart,'.dat']); %DAT file nameme for exporting stress

        % Define applied displacement in x direction
        displacementMagnitude = -0.5; %

        % linear material parameters
        E_youngs = 1666;
        v_poisson = 0.38;

        resultStruct.material.linear.E_youngs = E_youngs;
        resultStruct.material.linear.v_poisson = v_poisson;

        tolDir = 0.05;

        %(*@\codesubsection{Create TPMS}{lattice-create-TPMS}@*)
        [S,X,Y,Z] = gradTPMS(inputStruct);

        % Isosurface
        [F,V] = isosurface(X,Y,Z,S,levelset);

        %Capping ends
        [fc,vc]=isocaps(X,Y,Z,S,levelset, 'above');

        % Join, merge, and clean unused
        [f,v,c] = FV_arrange(F,V,fc,vc);

        resultStruct.f = f; resultStruct.v = v; resultStruct.c = c;
        %(*@\codesubsection{Visualise Surface}{lattice-visualise-surface}@*)
        cFigure;
        if figureTitles
            title('Face Labelling');
        end
        hp1 = gpatch(f,v,c,'none', 1);
        hp1.FaceColor = 'Flat';
        colormap(gca,gjet(6));
        axisGeom(gca,fontSize); camlight headlight;

        cFigure;
        if figureTitles
            title('Isosurface Orthographic Projection');
        end
        % Create the patch with uniform color
        hp2 = gpatch(f, v, 'none', 1);
        hp2.FaceColor = [0 0.4470 0.7410];
        hp2.EdgeColor = [0 0.4470 0.7410];
        axis equal;
        xlim([min(v(:,1)), max(v(:,1))]) % Adjust x-axis limits to your data range
        ylim([min(v(:,2)), max(v(:,2))]) % Adjust y-axis limits to your data range

        cFigure;
        if figureTitles
            title('Iso-surface');
        end
        hp3 = gpatch(f,v,c,'none', 1);
        hp3.FaceColor = [0 0.4470 0.7410];
        colormap(gca,gjet(6));
        axisGeom(gca,fontSize); camlight headlight;

        %(*@\codesubsection{Remesh Using geomgram}{lattice-remesh-using-geomgram}@*)
        optionStruct.pointSpacing=pointSpacing;
        [F,V]=ggremesh(f,v,optionStruct);
        C=zeros(size(F,1),1);

        resultStruct.F = F; resultStruct.V = V;

        % Visualizing geometry
        cFigure; hold on;

        if figureTitles
            title('Geogram remeshed','FontSize',fontSize);
        end

        gpatch(F,V,'w','k',1);
        axisGeom(gca,fontSize);
        camlight headlight;
        drawnow;

        %(*@\codesubsection{TetGen Tetrahedral Meshing}{lattice-tetgen-tetrahedral-meshing}@*)
        maxIterations = 20;
        iterationCount = 0;

        while iterationCount < maxIterations

            if coarseMesh
                optionStruct.pointSpacing = 0.4;
            end

            [F,V]=ggremesh(f,v,optionStruct);
            C=zeros(size(F,1),1);

            inputStruct.stringOpt='-pq1.2AaY';
            inputStruct.Faces=F;
            inputStruct.Nodes=V;
            inputStruct.holePoints=[];
            inputStruct.faceBoundaryMarker=C; %Face boundary markers
            inputStruct.regionPoints=getInnerPoint(F,V); %region points
            inputStruct.regionA=2*tetVolMeanEst(F,V);
            inputStruct.minRegionMarker=2; %Minimum region marker

            [meshOutput]=runTetGen(inputStruct); %Run tetGen

            numElements = size(meshOutput.elements, 1);

            fprintf("%i Elements generated\n",numElements);

            if coarseMesh
                break;
            elseif numElements > 1000000 % increase pointSpacing
                optionStruct.pointSpacing = optionStruct.pointSpacing + 0.1
            elseif numElements > 500000
                optionStruct.pointSpacing = optionStruct.pointSpacing + 0.05
            elseif numElements > 300000
                optionStruct.pointSpacing = optionStruct.pointSpacing + 0.03
            elseif numElements > 250000
                optionStruct.pointSpacing = optionStruct.pointSpacing + 0.01
            elseif numElements < 200000 % reduce pointSpacing
                optionStruct.pointSpacing = optionStruct.pointSpacing - 0.01
            else
                break;
            end

            iterationCount = iterationCount + 1;
        end

        if iterationCount == maxIterations
            warning('Maximum iterations reached without achieving desired element count.');
        end

        resultStruct.numElements = numElements;
        resultStruct.meshOutput = meshOutput;

        % Access model element and patch data
        Fb=meshOutput.facesBoundary;
        Cb=meshOutput.boundaryMarker;
        V=meshOutput.nodes;
        CE=meshOutput.elementMaterialID;
        E=meshOutput.elements;

        V_mesh = sum(tetVol(E,V,0)); % Mesh volume (through sum of tet elements volume)
        Infill_percentage = (V_mesh / (x_length * y_length * z_length) * 100); % assuming volume of 3

        %(*@\codesubsection{Visualise Mesh}{lattice-visualise-mesh}@*)
        meshView(meshOutput);

        %(*@\codesubsection{Node Labels and Selection Logic}{lattice-node-labels-and-selection-logic}@*)
        C_vertex=zeros(size(V,1),1);

        logic1=V(:,1)>(x_length-tolDir);
        logic2=V(:,1)<tolDir;

        C_vertex(logic1)=max(C_vertex(:))+1;
        C_vertex(logic2)=max(C_vertex(:))+1;

        logic3=V(:,3)>(z_length-tolDir)&V(:,1)>tolDir;
        logic4=V(:,3)<tolDir&V(:,1)>tolDir;

        C_vertex(logic3)=max(C_vertex(:))+1;
        C_vertex(logic4)=max(C_vertex(:))+1;

        logic5 = (V(:,2) > (y_length - tolDir)) & (V(:,1) > tolDir);
        logic6 = (V(:,2) < tolDir) & (V(:,1) > tolDir);

        C_vertex(logic5) = max(C_vertex(:)) + 1;
        C_vertex(logic6) = max(C_vertex(:)) + 1;

        resultStruct.C_vertex = C_vertex;

        % Visualizing vertex/node labels
        hf=cFigure;

        if figureTitles
            title('Boundary Nodes','FontSize',fontSize);
        end

        xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
        hold on;

        gpatch(Fb,V,'w','none',1);
        scatterV(V,10,C_vertex,'filled');

        axisGeom(gca,fontSize);
        colormap gjet; icolorbar;
        camlight headlight;
        drawnow;

        %(*@\codesubsection{Boundary Conditions}{lattice-boundary-conditions}@*)
        bcEncastreList=find(C_vertex==1); % Encastre vertices
        bcLoadList=find(C_vertex==2); % Load vertices
        bcTrackZ1=find(C_vertex==3);
        bcTrackZ0=find(C_vertex==4);
        bcTrackY0=find(C_vertex==5);
        bcTrackY1=find(C_vertex==6);

        bcTrackZ = union(bcTrackZ1, bcTrackZ0);
        bcTrackY = union(bcTrackY1, bcTrackY0);

        resultStruct.bcSets.bcEncastreList = bcEncastreList;
        resultStruct.bcSets.bcLoadList = bcLoadList;
        resultStruct.bcSets.bcTrackY = bcTrackY;
        resultStruct.bcSets.bcTrackY0 = bcTrackY0;
        resultStruct.bcSets.bcTrackY1 = bcTrackY1;

        resultStruct.bcSets.bcTrackZ = bcTrackZ;
        resultStruct.bcSets.bcTrackZ0 = bcTrackZ0;
        resultStruct.bcSets.bcTrackZ1 = bcTrackZ1;

        % Visualizing boundary conditions.
        hl=cFigure;
        if figureTitles
            title('Boundary Conditions','FontSize',fontSize);
        end
        xlabel('X','FontSize',fontSize); ylabel('Y','FontSize',fontSize); zlabel('Z','FontSize',fontSize);
        hold on;
        gpatch(Fb,V,'kw','none',0.5);

        hl(1)=plotV(V(bcEncastreList,:),'r*','MarkerSize',markerSize);
        hl(2)=plotV(V(bcLoadList,:),'k.','MarkerSize',markerSize);
        hl(3)=plotV(V(bcTrackZ,:),'.','MarkerSize',markerSize/2, 'MarkerEdgeColor', blue);
        hl(4)=plotV(V(bcTrackY,:),'c.','MarkerSize',markerSize/2, 'MarkerEdgeColor', cyan);
        legend(hl,{'BC Load','BC Encastre','BC Track Z','BC Track Y'});

        axisGeom(gca,fontSize);
        camlight headlight;
        drawnow;

        %(*@\codesubsection{Save All Figures}{lattice-save-all-figures}@*)
        figsavePath = fullfile(savePath, 'Figures');
        mkdir(figsavePath);
        nFigures = saveFigures(savePath,doSaveFigures,nFigures);

        %(*@\codesubsection{Abaqus Input Structure}{lattice-abaqus-input-structure}@*)
        %%--> Heading
        abaqus_spec.Heading.COMMENT{1}='Job name: ABAQUS inp file creation demo';
        abaqus_spec.Heading.COMMENT{2}='Generated by: GIBBON';

        %%--> Preprint
        abaqus_spec.Preprint.ATTR.echo='NO';
        abaqus_spec.Preprint.ATTR.model='NO';
        abaqus_spec.Preprint.ATTR.history='NO';
        abaqus_spec.Preprint.ATTR.contact='NO';
        abaqus_spec.Preprint.ATTR.contact='NO';

        %--> Part

        % Node
        nodeIds=(1:1:size(V,1))';
        abaqus_spec.Part.COMMENT='This section defines the part geometry in terms of nodes and elements';
        abaqus_spec.Part.ATTR.name='Lattice';
        abaqus_spec.Part.Node={nodeIds,V};

        % Element
        elementIds=(1:1:size(E,1));
        abaqus_spec.Part.Element{1}.ATTR.type='C3D4';
        abaqus_spec.Part.Element{1}.VAL={elementIds(:),E};

        % Element sets
        abaqus_spec.Part.Elset{1}.ATTR.elset='Set-1';
        abaqus_spec.Part.Elset{1}.VAL=elementIds;

        abaqus_spec.Part.Solid_section.ATTR.elset='Set-1';
        abaqus_spec.Part.Solid_section.ATTR.material='Elastic';

        %%--> Assembly
        abaqus_spec.Assembly.ATTR.name='Assembly-1';
        abaqus_spec.Assembly.Instance.ATTR.name='Lattice-assembly';
        abaqus_spec.Assembly.Instance.ATTR.part='Lattice';

        abaqus_spec.Assembly.Nset{1}.ATTR.nset='Set-1';
        abaqus_spec.Assembly.Nset{1}.ATTR.instance='Lattice-assembly';
        abaqus_spec.Assembly.Nset{1}.VAL=bcEncastreList(:)';

        abaqus_spec.Assembly.Nset{2}.ATTR.nset='Set-2';
        abaqus_spec.Assembly.Nset{2}.ATTR.instance='Lattice-assembly';
        abaqus_spec.Assembly.Nset{2}.VAL=bcLoadList(:)';

        abaqus_spec.Assembly.Nset{3}.ATTR.nset='all';
        abaqus_spec.Assembly.Nset{3}.ATTR.instance='Lattice-assembly';
        abaqus_spec.Assembly.Nset{3}.VAL=1:1:size(V,1);

        %%--> Material
        abaqus_spec.Material.ATTR.name='Elastic';
        abaqus_spec.Material.Elastic=[E_youngs v_poisson];

        %%--> Step
        abaqus_spec.Step.ATTR.name='Step-1';
        abaqus_spec.Step.ATTR.nlgeom='YES';
        abaqus_spec.Step.Static=[0.01 1 1e-5 0.1];

        %--> Boundary
        abaqus_spec.Step.Boundary{1}.VAL={'Set-1', 'XSYMM'};
        abaqus_spec.Step.Boundary{2}.VAL={'Set-2',[1,1],displacementMagnitude};

        %--> Output
        %Nodal coordinates
        abaqus_spec.Step.Node_print{1}.ATTR.nset='all';
        abaqus_spec.Step.Node_print{1}.ATTR.frequency = 1;
        abaqus_spec.Step.Node_print{1}.VAL='COORD';

        abaqus_spec.Step.El_print{1}.VAL='S';
        abaqus_spec.Step.El_print{2}.VAL='E';

        %(*@\codesubsection{Writing INP File}{lattice-writing-INP-file}@*)
        if updateSim

            abaqusStruct2inp(abaqus_spec,abaqusInpFileName);

            %(*@\codesubsection{Run Abaqus Job}{lattice-run-abaqus-job}@*)
            lockFileName=fullfile(savePath,[abaqusInpFileNamePart,'.lck']);

            if exist(lockFileName,'file')
                warning('Lockfile found and deleted')
                delete(lockFileName);
            end

            oldPath=pwd; %Get current working directory
            cd(savePath); %Set new working directory to match save path

            abaqusPath='C:\SIMULIA\Commands\abaqus.bat';%'/usr/bin/abaqus'; %Abaqus excute command or path
            [runFlag, cmdOut] = system([abaqusPath,' inp=',abaqusInpFileName,' job=',abaqusInpFileNamePart,' interactive ask_delete=OFF ', 'cpus=',int2str(4)]);
            disp(cmdOut);

            cd(oldPath); %Restore working directory
        else
            fprintf('Abaqus files for [%.1g] Lattice density exist and should not be overwritten, proceeding\n', i);
        end

        %(*@\codesubsection{Import and Visualise Abaqus Results}{lattice-import-and-visualise-abaqus-results}@*)
        saveFile = fullfile(savePath, 'simulation_results.mat');

        try
            % Attempt to import the Abaqus data file.
            if loadSimDataFromMat && isfile(saveFile)
                try
                    load(saveFile);
                    abaqusData = resultStruct.abaqusData;
                catch
                    abaqusData = importAbaqusDat(abaqusDATFileName);
                end
            else
                abaqusData = importAbaqusDat(abaqusDATFileName);
            end
        catch ME
            rethrow(ME)
            error('Check Abaqus License, connect to campus wifi or campus VPN')
        end

        %(*@\codesubsection{Fetch Element Data}{lattice-fetch-element-data}@*)
        E_effectiveStress=zeros(size(E,1),numel(abaqusData.STEP.INCREMENT)+1);
        E_effectiveStrain=zeros(size(E,1),numel(abaqusData.STEP.INCREMENT)+1);
        for q=1:1:numel(abaqusData.STEP.INCREMENT)

            for d=1:1:2
                dataSet=abaqusData.STEP.INCREMENT(q).elementOutput(d); %Get data structure
                switch d
                    case 1
                        % Get element data across all 1 integration points
                        D11_PT=reshape(dataSet.data.S11,1,size(E,1))';
                        D22_PT=reshape(dataSet.data.S22,1,size(E,1))';
                        D33_PT=reshape(dataSet.data.S33,1,size(E,1))';
                        D12_PT=reshape(dataSet.data.S12,1,size(E,1))';
                        D13_PT=reshape(dataSet.data.S13,1,size(E,1))';
                        D23_PT=reshape(dataSet.data.S23,1,size(E,1))';
                    case 2
                        % Get element data across all 1 integration points
                        D11_PT=reshape(dataSet.data.E11,1,size(E,1))';
                        D22_PT=reshape(dataSet.data.E22,1,size(E,1))';
                        D33_PT=reshape(dataSet.data.E33,1,size(E,1))';
                        D12_PT=reshape(dataSet.data.E12,1,size(E,1))';
                        D13_PT=reshape(dataSet.data.E13,1,size(E,1))';
                        D23_PT=reshape(dataSet.data.E23,1,size(E,1))';
                end

                % Get mean element metrics
                D11=mean(D11_PT,2);
                D22=mean(D22_PT,2);
                D33=mean(D33_PT,2);
                D12=mean(D12_PT,2);
                D13=mean(D13_PT,2);
                D23=mean(D23_PT,2);

                % Calculate effective (Von Mises) metrics
                switch d
                    case 1
                        E_effectiveStress(:,q+1)=(1/sqrt(2)).*sqrt((D11-D22).^2+(D11-D33).^2+(D33-D11).^2+6*D12.^2+6*D23.^2+6*D13.^2);
                    case 2
                        E_effectiveStrain(:,q+1)=(1/sqrt(2)).*sqrt((D11-D22).^2+(D11-D33).^2+(D33-D11).^2+6*D12.^2+6*D23.^2+6*D13.^2);
                end
            end
        end

        %(*@\codesubsection{Animated Deformations}{lattice-animated-deformations}@*)
        V_def=[abaqusData.STEP(1).INCREMENT(end).nodeOutput(1).data.COOR1...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput(1).data.COOR2...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput(1).data.COOR3];
        U=V_def-V; %Displacements

        colorDataVertices=sqrt(sum(U.^2,2)); %Displacement magnitude data for coloring

        timeVec=[0 abaqusData.STEP(1).INCREMENT(:).TOTAL_TIME_COMPLETED];

        % Get limits for plotting
        minV=min([V;V_def],[],1); %Minima
        maxV=max([V;V_def],[],1); %Maxima

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure
        gtitle([abaqusInpFileNamePart,': Displacement data']);
        hp=gpatch(Fb,V_def,colorDataVertices,'none',1); %Add graphics object to animate
        gpatch(Fb,V,0.5*ones(1,3),'k',0.25); %A static graphics object
        axisGeom(gca,fontSize);
        colormap(gjet(250)); colorbar;
        caxis([0 max(colorDataVertices)]);
        axis([minV(1) maxV(1) minV(2) maxV(2) minV(3) maxV(3)]); %Set axis limits statically
        view(130,25); %Set view direction
        camlight headlight;

        % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:numel(timeVec) %Loop over time increments
            if qt>1
                V_def=[abaqusData.STEP(1).INCREMENT(qt-1).nodeOutput.data.COOR1...
                    abaqusData.STEP(1).INCREMENT(qt-1).nodeOutput.data.COOR2...
                    abaqusData.STEP(1).INCREMENT(qt-1).nodeOutput.data.COOR3];
            else
                V_def=V;
            end
            U=V_def-V; %Displacements
            colorDataVertices=sqrt(sum(U(:,3).^2,2)); %New color data

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_def,colorDataVertices}; %Property values for to set in order to animate
        end

        anim8(hf,animStruct); %Initiate animation feature
        drawnow;

        %(*@\codesubsection{Animated Von Mises}{lattice-animated-von-mises}@*)
        V_def=[abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR1...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR2...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR3];
        U=V_def-V; %Displacements

        [FE,C_FE_effectiveStress]=element2patch(E,E_effectiveStress(:,end),'tet4');
        [indBoundary]=tesBoundary(FE,V);
        FEb=FE(indBoundary,:);
        colorDataFaces=C_FE_effectiveStress(indBoundary);

        timeVec=[0 abaqusData.STEP(1).INCREMENT(:).TOTAL_TIME_COMPLETED];

        % Get limits for plotting
        minV=min([V;V_def],[],1); %Minima
        maxV=max([V;V_def],[],1); %Maxima

        % Create basic view and store graphics handle to initiate animation
        hf=cFigure; %Open figure
        gtitle([abaqusInpFileNamePart,': Effective Von Mises Stress (MPa)']);
        hp=gpatch(FEb,V_def,colorDataFaces,'none',1); %Add graphics object to animate
        gpatch(FEb,V,0.5*ones(1,3),'k',0.25); %A static graphics object
        axisGeom(gca,fontSize);
        colormap(gjet(250)); colorbar;
        caxis([0 max(E_effectiveStress(:))]);
        axis([minV(1) maxV(1) minV(2) maxV(2) minV(3) maxV(3)]); %Set axis limits statically
        view(130,25); %Set view direction
        camlight headlight;

        % Set up animation features
        animStruct.Time=timeVec; %The time vector
        for qt=1:1:numel(timeVec) %Loop over time increments
            if qt>1
                V_def=[abaqusData.STEP(1).INCREMENT(qt-1).nodeOutput.data.COOR1...
                    abaqusData.STEP(1).INCREMENT(qt-1).nodeOutput.data.COOR2...
                    abaqusData.STEP(1).INCREMENT(qt-1).nodeOutput.data.COOR3];
            else
                V_def=V;
            end
            U=V_def-V; %Displacements

            [~,C_FE_effectiveStress]=element2patch(E,E_effectiveStress(:,qt),'tet4');
            colorDataFaces=C_FE_effectiveStress(indBoundary); %New color data

            %Set entries in animation structure
            animStruct.Handles{qt}=[hp hp]; %Handles of objects to animate
            animStruct.Props{qt}={'Vertices','CData'}; %Properties of objects to animate
            animStruct.Set{qt}={V_def,colorDataFaces}; %Property values for to set in order to animate
        end

        anim8(hf,animStruct); %Initiate animation feature
        drawnow;

        %(*@\codesubsection{Count Distorted Elements}{count-distorted-elements}@*)
        fid = fopen(abaqusDATFileName, 'r');
        if fid == -1
            error('Unable to open the file: %s', abaqusDATFileName);
        end

        distortedElements = [];
        pattern = 'WARNING:\s*(\d+)\s*elements';  % regex pattern to capture the number

        % Read file line-by-line until we find the matching warning string
        while ~feof(fid)
            tline = fgetl(fid);
            if contains(tline, '***WARNING:')
                tokens = regexp(tline, pattern, 'tokens');
                if ~isempty(tokens)
                    distortedElements = str2double(tokens{1}{1});
                    break;
                end
            end
        end

        fclose(fid);

        if isempty(distortedElements)
            warning('The warning string was not found in the file.');
        end

        resultStruct.distortedElements = distortedElements;

        %(*@\codesubsection{Poisson's Calculated at All Increments}{lattice-poissons-calculated-at-all-increments}@*)
        numSteps = size(abaqusData.STEP(1).INCREMENT,2); % Get number of increments

        % Preallocate arrays with NaNs
        strain_z_mean = zeros(numSteps, 1); strain_z_median = zeros(numSteps, 1); strain_z_max = zeros(numSteps, 1);
        strain_y_mean = zeros(numSteps, 1); strain_y_median = zeros(numSteps, 1); strain_y_max = zeros(numSteps, 1);

        strain_x = zeros(numSteps, 1); U_x = zeros(numSteps, 1);

        for i = 1:numSteps

            %Getting nodal coordinates
            V_def=[abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR1...
                abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR2...
                abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR3];
            U=V_def-V; %Displacements

            displacementMagnitude = U(bcLoadList(1),1); % checks how much single node was displaced in x in this timestep
            U_x(i) = displacementMagnitude;
            strain_x(i) = ((-displacementMagnitude + x_length) - x_length) / x_length; % = deformed length / initial length

            strain_z_mean(i) = (mean(U(bcTrackZ1,3)) - mean(U(bcTrackZ0,3))) / z_length;
            strain_z_median(i) = (median(U(bcTrackZ1,3)) - median(U(bcTrackZ0,3))) / z_length;
            strain_z_max(i) = (min(U(bcTrackZ1,3)) - max(U(bcTrackZ0,3))) / z_length;

            % Y displacement of tracked sides
            strain_y_mean(i) = (mean(U(bcTrackY0,2)) - mean(U(bcTrackY1,2))) / y_length;
            strain_y_median(i) = (median(U(bcTrackY0,2)) - median(U(bcTrackY1,2))) / y_length;
            strain_y_max(i) = (min(U(bcTrackY0,2)) - max(U(bcTrackY1,2))) / y_length;
        end

        poisson_xz_mean = -strain_z_mean ./ strain_x; % note that since strain_y = U_y_mean / 1 we can conclude for this case strain_y = U_y_mean
        poisson_xz_median = -strain_z_median ./ strain_x;
        poisson_xz_max = -strain_z_max ./ strain_x;

        poisson_xy_mean = -strain_y_mean ./ strain_x;
        poisson_xy_median = -strain_y_median ./ strain_x;
        poisson_xy_max = -strain_y_max ./ strain_x;

        %(*@\codesubsection{Z Displacement with Step Increments}{lattice-z-displacement-with-step-increments}@*)
        nRows = ceil(sqrt(numSteps));
        nCols = ceil(numSteps / nRows);

        cFigure;
        if figureTitles
            gtitle('Tracked Z Displacement (mm)')
        end
        for i = numSteps:-1:1

            h = subplot(nRows, nCols, i);

            V_def=[abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR1...
                abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR2...
                abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR3];
            U = V_def-V; %Displacements

            displacementMagnitude = U(bcLoadList(1),1);

            xz0 = abs(V(bcTrackZ0,1));
            Uz0 = abs(U(bcTrackZ0,3));
            xz1 = abs(V(bcTrackZ1,1));
            Uz1 = abs(U(bcTrackZ1,3));

            if i == numSteps
                globalYMax = max(vertcat(Uz0,Uz1));
            end

            % Plot the original data points
            hold on;
            scatter(xz0, Uz0, 'x', 'MarkerEdgeColor', blue);
            scatter(xz1, Uz1, 'x', 'MarkerEdgeColor', cyan);

            % Labels and legend
            xlabel('Global X Coordinate (mm)');
            ylabel('Abs Z Displacement (mm)');
            hTitle = title(sprintf('%.3g X (mm) - Step %i/%i', displacementMagnitude,i,numSteps));
            set(hTitle, 'FontSize', 8);
            grid on;

            ylim([0, globalYMax]);
            xlim([0, x_length]);

        end
        drawnow

        %(*@\codesubsection{Y Displacement with Step Increments}{lattice-y-displacement-with-step-increments}@*)
        cFigure;
        if figureTitles
            gtitle('Tracked Y Displacement (mm)')
        end
        for i = numSteps:-1:1

            h = subplot(nRows, nCols, i);

            % subplot % SUBPLOT with numstep thingy in here
            V_def=[abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR1...
                abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR2...
                abaqusData.STEP(1).INCREMENT(i).nodeOutput.data.COOR3];
            U = V_def-V; %Displacements

            displacementMagnitude = U(bcLoadList(1),1);

            xy0 = abs(V(bcTrackY0,1));
            Uz0 = abs(U(bcTrackY0,2));
            xy1 = abs(V(bcTrackY1,1));
            Uz1 = abs(U(bcTrackY1,2));

            if i == numSteps
                globalYMax = max(vertcat(Uz0,Uz1));
            end

            % Plot the original data points
            hold on;
            scatter(xy0, Uz0, 'x', 'MarkerEdgeColor', blue);
            scatter(xy1, Uz1, 'x', 'MarkerEdgeColor', cyan);

            % Labels and legend
            xlabel('Global X Coordinate (mm)');
            ylabel('Abs Y Displacement (mm)');
            hTitle = title(sprintf('%.3g X (mm) - Step %i/%i', displacementMagnitude,i,numSteps));
            set(hTitle, 'FontSize', 8);
            grid on;

            ylim([0, globalYMax]);
            xlim([0, x_length]);

        end
        drawnow;

        %(*@\codesubsection{Final Step Tracked Z Displacement}{lattice-final-step-tracked-z-displacement}@*)
        cFigure;
        if figureTitles
            gtitle('Tracked Z Displacement (mm)')
        end

        hold on;

        V_def=[abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR1...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR2...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR3];
        U = V_def-V; %Displacements

        displacementMagnitude = U(bcLoadList(1),1);

        xz0 = abs(V(bcTrackZ0,1));
        Uz0 = abs(U(bcTrackZ0,3));
        xz1 = abs(V(bcTrackZ1,1));
        Uz1 = abs(U(bcTrackZ1,3));

        resultStruct.xz0 = xz0; resultStruct.Uz0 = Uz0;
        resultStruct.xz1 = xz1; resultStruct.Uz1 = Uz1;
        scatter(xz0, Uz0, 'x', 'MarkerEdgeColor', blue, 'LineWidth', 1.5);
        scatter(xz1, Uz1, 'x', 'MarkerEdgeColor', cyan, 'LineWidth', 1.5);

        xlabel('Global X Coordinate (mm)');
        ylabel('Absolute Z Displacement (mm)');
        grid on; hold off;

        %(*@\codesubsection{Final Step Tracked Y Displacement}{lattice-final-step-tracked-y-displacement}@*)
        cFigure;
        if figureTitles
            gtitle('Tracked Y Displacement (mm)')
        end

        hold on;

        V_def=[abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR1...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR2...
            abaqusData.STEP(1).INCREMENT(end).nodeOutput.data.COOR3];
        U = V_def-V; %Displacements

        xy0 = abs(V(bcTrackY0,1));
        Uy0 = abs(U(bcTrackY0,2));
        xy1 = abs(V(bcTrackY1,1));
        Uy1 = abs(U(bcTrackY1,2));

        resultStruct.xy0 = xy0; resultStruct.Uy0 = Uy0;
        resultStruct.xy1 = xy1; resultStruct.Uy1 = Uy1;
        scatter(xy0, Uy0, 'x', 'MarkerEdgeColor', blue, 'LineWidth', 1.5);
        scatter(xy1, Uy1, 'x', 'MarkerEdgeColor', cyan, 'LineWidth', 1.5);

        xlabel('Global X Coordinate (mm)');
        ylabel('Absolute Y Displacement (mm)');
        grid on;

        %(*@\codesubsection{Poisson's Based on Displacement Graph}{lattice-poisson-based-on-displacement-graph}@*)
        resultStruct.U_x = U_x;
        resultStruct.poisson.poisson_xz_mean = poisson_xz_mean;
        resultStruct.poisson.poisson_xz_median = poisson_xz_median;
        resultStruct.poisson.poisson_xy_mean = poisson_xy_mean;
        resultStruct.poisson.poisson_xy_median = poisson_xy_median;

        cFigure;
        hold on; grid on;

        if figureTitles
            gtitle('Poisson''s Ratio throughout all simulation steps')
        end

        % Plot poisson in xz (blue)
        plot(U_x, poisson_xz_mean, '-', 'LineWidth', 2,'Color',blue);
        plot(U_x, poisson_xz_median, '--', 'LineWidth', 2,'Color',blue);

        % Plot poisson in xy (cyan)
        plot(U_x, poisson_xy_mean, '-', 'LineWidth', 2,'Color',cyan);
        plot(U_x, poisson_xy_median, '--', 'LineWidth', 2,'Color',cyan);

        ylabel("Poisson's Ratio", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal")
        xlabel("Applied X Displacement (mm)", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal")

        legend({'poisson\_xz\_mean', 'poisson\_xz\_median', ...
            'poisson\_xy\_mean', 'poisson\_xy\_median'}, ...
            'Location', 'northeast');

        set(gca, 'XDir', 'reverse');
        % Final poisson value (average from all step increments)

        %(*@\codesubsection{Final Poisson Value}{lattice-final-poisson-value}@*)
        % Final poisson value taken as an average from all step increments
        % In reality previous plot did show poisson does change very slightly with large displacements
        poisson_xz_mean = mean(poisson_xz_mean);
        poisson_xz_median = mean(poisson_xz_median);

        poisson_xy_mean = mean(poisson_xy_mean);
        poisson_xy_median = mean(poisson_xy_median);

        resultStruct.poisson_xz_mean = poisson_xz_mean;
        resultStruct.poisson_xz_median = poisson_xz_median;
        resultStruct.poisson_xy_mean = poisson_xy_mean;
        resultStruct.poisson_xy_median = poisson_xy_median;

        appliedForce = zeros(numSteps, 1);
        stress = zeros(numSteps, 1);
        encastreTotalArea = zeros(numSteps, 1);

        for i = 1:numSteps
            stress(i) = mean(E_effectiveStress(bcEncastreList,i));

            encastreFaces = F(bcEncastreList, :);
            encastreFaceAreas = patchArea(encastreFaces, V);
            encastreTotalArea(i) = sum(encastreFaceAreas);

        end

        appliedForce = stress .* encastreTotalArea; % N

        resultStruct.appliedForce = appliedForce;
        resultStruct.encastreTotalArea = encastreTotalArea;

        cFigure;
        hold on; grid on;

        if figureTitles
            gtitle('Young''s Modulus')
        end

        resultStruct.strain.strain_y_mean = strain_y_mean;
        resultStruct.strain.strain_y_median = strain_y_median;
        resultStruct.strain.strain_z_mean = strain_z_mean;
        resultStruct.strain.strain_z_median = strain_z_median;
        resultStruct.stress = stress;

        plot(abs(strain_y_mean),stress,'-','Color',blue);
        plot(abs(strain_y_median),stress,'--','Color',blue);
        plot(abs(strain_z_mean),stress,'-','Color',cyan);
        plot(abs(strain_z_median),stress,'--','Color',cyan);

        legend('Mean Young''s Modulus xy', 'Median Young''s Modulus xy', 'Mean Young''s Modulus xz', 'Median Young''s Modulus xz','Location','best');

        xlabel("Strain", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal")
        ylabel("Stress (MPa)", "FontName", "Helvetica", "FontAngle", "normal", "FontWeight", "normal")

        %(*@\codesubsection{Young's Moduli}{lattice-youngs-moduli}@*)
        youngsModulus = zeros(1,4);
        R_squared = zeros(1,4);

        for i = 1:4

            if i == 1
                x = abs(strain_z_mean);
            elseif i == 2
                x = abs(strain_z_median);
            elseif i == 3
                x = abs(strain_y_mean);
            elseif i == 4
                x = abs(strain_y_median);
            end

            y = stress;

            [coefficients, stats] = polyfit(x,y, 1);
            R_squared(i) = stats.rsquared;

            youngsModulus(i) = coefficients(1);
        end

        resultStruct.youngs_modulus.youngs_modulus_xz_mean = youngsModulus(1);
        resultStruct.youngs_modulus.youngs_modulus_xz_median = youngsModulus(2);
        resultStruct.youngs_modulus.youngs_modulus_xy_mean = youngsModulus(3);
        resultStruct.youngs_modulus.youngs_modulus_xy_median = youngsModulus(4);

        resultStruct.youngs_modulus.youngs_modulus_xz_mean_R_squared = R_squared(1);
        resultStruct.youngs_modulus.youngs_modulus_xz_median_R_squared = R_squared(2);
        resultStruct.youngs_modulus.youngs_modulus_xy_mean_R_squared = R_squared(3);
        resultStruct.youngs_modulus.youngs_modulus_xy_median_R_squared = R_squared(4);

        %(*@\codesubsection{Re-Save New Figures}{lattice-resave-new-figures}@*)
        nFigures = saveFigures(savePath,doSaveFigures,nFigures);

        %(*@\codesubsection{Save to Struct}{lattice-save-to-struct}@*)
        resultStruct.maxVonMises = max(E_effectiveStress(:,end));
        resultStruct.percentile95VonMises = prctile(E_effectiveStress(:,end), 95);

        resultStruct.abaqusData = abaqusData;
        resultStruct.meshOutput = meshOutput;
        resultStruct.timestamp = datestr(now, 'yyyy-mm-dd HH:MM:SS'); % Save execution time

        resultStruct.infill_percentage = Infill_percentage;
        resultStruct.mesh_volume = V_mesh;
        resultStruct.level_set = levelset;

        resultStruct.lattice.x_length = x_length;
        resultStruct.lattice.y_length = y_length;
        resultStruct.lattice.z_length = z_length;

        resultStruct.E_effectiveStrain = E_effectiveStrain;
        resultStruct.E_effectiveStress = E_effectiveStress;

        resultStruct.coarseMesh = coarseMesh;
        resultStruct.figureTitles = figureTitles;

        resultStruct.pointSpacing = optionStruct.pointSpacing;

        saveFile = fullfile(savePath, 'simulation_results.mat');

        % Save struct
        save(saveFile, 'resultStruct');

    end

    %(*@\codesubsection{Stop the Data Collector Set}{lattice-stop-the-data-collector-set}@*)
    system('logman stop CPU_Log');

catch ME
    system('logman stop CPU_Log');
    rethrow(ME);
end