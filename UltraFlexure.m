classdef UltraFlexure < handle
    properties
        % PDE model
        model = createpde('structural','modal-solid')

        % structural properties
        E = 68.9E9
        nu = 0.33
        rho = 2700

        % modal analysis
        min_freq = -1 % Hz
        max_freq = 20E3 % Hz

        % modal analysis results
        RF
        modeID
        tmodalResults
        frames
    end

    methods
        function self = UltraFlexure(model_name)
            importGeometry(self.model,model_name);
            structuralBC(self.model,'Face',[7,29],'Constraint','roller')
            structuralProperties(self.model,'YoungsModulus',self.E, ...
                           'PoissonsRatio',self.nu, ...
                           'MassDensity',self.rho);
            generateMesh(self.model,'Hmax',0.003);
        end

        function solveModal(self)
            self.RF = solve(self.model,'FrequencyRange',[self.min_freq,self.max_freq]*2*pi);
            self.modeID = 1:numel(self.RF.NaturalFrequencies);
        end

        function showModes(self)
            self.tmodalResults = table(self.modeID.',self.RF.NaturalFrequencies/2/pi);
            self.tmodalResults.Properties.VariableNames = {'Mode','Frequency'};
            disp(self.tmodalResults); 
        end

        function generateFrames(self, modes)
            % Animation parameters
            scale = 0.0005;
            nFrames = 30;
            flexibleModes = modes;
             
            % Create a model for plotting purpose.
            deformedModel = createpde('structural','modal-solid');
             
            % Undeformed mesh data
            nodes = self.RF.Mesh.Nodes;
            elements = self.RF.Mesh.Elements;
             
            % Construct pseudo time-vector that spans one period of first six flexible
            % modes.
            omega = self.RF.NaturalFrequencies(flexibleModes);
            timePeriod = 2*pi./self.RF.NaturalFrequencies(flexibleModes);
             
            h = figure('units','normalized','outerposition',[0 0 1 1]);
            % Plot deformed shape of the first six flexible modes and capture frame for
            % each pseudo time step.
            for n = 1:nFrames
                for modeID = 1:6
                    % Construct a modal deformation matrix and its magnitude.
                    modalDeformation = [self.RF.ModeShapes.ux(:,flexibleModes(modeID))';
                        self.RF.ModeShapes.uy(:,flexibleModes(modeID))';
                        self.RF.ModeShapes.uz(:,flexibleModes(modeID))'];
                    
                    modalDeformationMag = sqrt(modalDeformation(1,:).^2 + ...
                        modalDeformation(2,:).^2 + ...
                        modalDeformation(3,:).^2);
                    
                    % Compute nodal locations of deformed mesh.
                    pseudoTimeVector = linspace(0,timePeriod(modeID),nFrames);
                    nodesDeformed = nodes + scale.*modalDeformation*sin(omega(modeID).*pseudoTimeVector(n));
                    
                    % Construct a deformed geometric shape using displaced nodes and
                    % elements from unreformed mesh data.
                    geometryFromMesh(deformedModel,nodesDeformed,elements);
                    
                    % Plot the deformed mesh with magnitude of mesh as color plot.
                    subplot(2,3,modeID)
                    pdeplot3D(deformedModel,'ColorMapData', modalDeformationMag)
                    title(sprintf(['Flexible Mode %d\n', ...
                        'Frequency = %g Hz'], ...
                        modeID,omega(modeID)/2/pi));
                    
                    
                    % Remove axes triad and colorbar for clarity
                    colorbar off
                    axis([0.01,0.1,0.01,0.06,0.002,0.003]);
                    delete(findall(gca,'type','quiver'));
                    qt = findall(gca,'type','text');
                    set(qt(1:3),'Visible','off')
                    
                    % Remove deformed geometry to reuse to model for next mode.
                    deformedModel.Geometry = [];
                end
                % Capture a frame of six deformed mode for time instant
                frames(n) = getframe(h);
            end
            self.frames = frames;
        end

        function showMovie(self)
            movie(figure('units','normalized','outerposition',[0 0 1 1]),self.frames,5,30)
        end

    end
end