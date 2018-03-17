classdef FeaDataStructure < handle
    properties
    % constants
    end_thk = 5*1e-3
    wall_gap = 2.4*1e-3
    mirror_width = 5*1e-3
    press_len = 2*1e-3
    press_width = 0.5*1e-3
    total_thk = 7.5e-3
    top_key_height = 3e-3
    support_poly_gap = 1e-3

    % variables
    pinch_len
    press_pos
    wall_thk
    wall_len
    angle_wall_len
    theta

    % dependents
    angle_len
    total_width
    total_len

    % pde
    model
    edges
    pde_c
    pde_e = 1300e6 %2310e6 % Pa (reported between 1000 and 2000 MPa) % 1300 MPa from other datasheet?
    pde_gnu = .408 % Poisson's ratio of the material: 0.408
    pde_f = [0 0]' % No body forces
    % with press support poly
    geo_names = char('body','corner1','corner2','corner3','corner4',...
                     'left_poly','right_poly','bottom_poly','top_poly','press_support_poly',...
                     'top_key','bottom_press','top_press')';
    geo_formula = '(body+bottom_press+top_press+top_key)-corner1-corner2-corner3-corner4-left_poly-right_poly-bottom_poly-top_poly+press_support_poly';
    % without press support poly
    % geo_names = char('body','corner1','corner2','corner3','corner4',...
    %                  'left_poly','right_poly','bottom_poly','top_poly',...
    %                  'top_key','bottom_press','top_press')';
    % geo_formula = '(body+bottom_press+top_press+top_key)-corner1-corner2-corner3-corner4-left_poly-right_poly-bottom_poly-top_poly';
    force_app
    surface_traction
    results
    results_p
    results_e
    results_t
    deflection_x
    deflection_y
    theta_buffer
    force_app_buffer
    force_pos_buffer
    deflection_x_buffer
    deflection_y_buffer
    end

    methods
        function self = FeaDataStructure(wall_thk, wall_len, theta, force_app, pinch_len, press_pos)
            G = self.pde_e/(2.*(1+self.pde_gnu));
            % mu = 2*G*self.pde_gnu/(1-self.pde_gnu); % plane stress
            mu = 2*G*self.pde_gnu/(1-2*self.pde_gnu); % plane strain
            self.pde_c = [2*G+mu; 0; G;   0; G; mu; 0;  G; 0; 2*G+mu];
            self.setAndSolve(wall_thk, wall_len, theta, force_app, pinch_len, press_pos);
        end

        function setAndSolve(self, wall_thk, wall_len, theta, force_app, pinch_len, press_pos)
            self.setGeometry(wall_thk, wall_len, theta, pinch_len, press_pos);
            self.setModel(force_app);
            self.solveModel();
        end

        function setGeometry(self, wall_thk, wall_len, theta, pinch_len, press_pos)
            len_ratio = 1.0;
            thk_ratio = 1.0;
            taper_ratio = 1.0;
            self.wall_thk = wall_thk*1e-3;
            self.wall_len = wall_len*1e-3;
            self.pinch_len = pinch_len*1e-3;
            self.press_pos = press_pos*1e-3;
            self.angle_wall_len = self.wall_len*len_ratio;
            self.theta = theta;

            % dependents
            self.angle_len = self.angle_wall_len*sin(deg2rad(theta));
            self.total_width = 4*self.wall_thk+2*self.wall_gap+2*self.angle_len+self.mirror_width;
            self.total_len = 2*self.end_thk+2*self.wall_len+self.pinch_len;

            % csg geometry for basic geometry
            main_rect = [3,4,0,0,self.total_len,self.total_len,0,self.total_width,self.total_width,0];
            corner_rect = [3,4,self.end_thk,self.end_thk,self.end_thk+self.wall_len,self.end_thk+self.wall_len,self.wall_thk,self.wall_thk+self.wall_gap,self.wall_thk+self.wall_gap,self.wall_thk]; % bottom left rect
            x_off = self.wall_len+self.pinch_len;
            y_off = 2*self.wall_thk+self.wall_gap+2*self.angle_len+self.mirror_width;
            corner_rects = [corner_rect; % all four rects in columns
                            corner_rect+[0,0,0,0,0,0,y_off,y_off,y_off,y_off];
                            corner_rect+[0,0,x_off,x_off,x_off,x_off,0,0,0,0];
                            corner_rect+[0,0,x_off,x_off,x_off,x_off,y_off,y_off,y_off,y_off]];
            x1 = self.end_thk;
            x2 = x1+self.wall_len;
            x3 = x2+self.angle_wall_len*cos(deg2rad(self.theta));
            x4 = self.end_thk+self.wall_len+self.pinch_len;
            x5 = x4+self.angle_wall_len*cos(deg2rad(self.theta));
            x6 = x4+self.wall_len;
            x7 = x2+self.wall_thk/sin(deg2rad(self.theta));%*taper_ratio;
            x8 = x3+self.wall_thk/sin(deg2rad(self.theta));
            x9 = x5-self.wall_thk/sin(deg2rad(self.theta))*thk_ratio;
            x10 = x4-self.wall_thk/sin(deg2rad(self.theta))*thk_ratio;%*taper_ratio;
            y1 = 2*self.wall_thk+self.wall_gap;
            y2 = y1+2*self.angle_len+self.mirror_width;
            y3 = y2-self.angle_len;
            y4 = y3-self.mirror_width;
            left_poly = [2,6,x1,x1,x2,x3,x3,x2,y1,y2,y2,y3,y4,y1];
            right_poly = [2,6,x4,x5,x5,x4,x6,x6,y1,y4,y3,y2,y2,y1];
            mid_poly = [2,4,x7,x8,x9,x10,y1,y4,y4,y1];
            poly_y_off = self.angle_len+self.mirror_width;
            mid_polys = [mid_poly;
                         mid_poly(1:6) [mid_poly(9:10) mid_poly(7:8)]+poly_y_off];
            press_x1 = .5*self.total_len-.5*self.pinch_len;
            press_x2 = press_x1+self.pinch_len;
            press_x3 = .5*self.total_len-.5*self.press_len+self.press_pos;
            press_x4 = press_x3+self.press_len;
            press_y1 = 0;
            press_y2 = -self.press_width;
            top_key_y1 = self.total_width;
            top_key_y2 = self.total_width+self.top_key_height;
            press_y3 = top_key_y2;
            press_y4 = press_y3+self.press_width;
            bottom_press_rect = [3,4,press_x1,press_x1,press_x2,press_x2,press_y1,press_y2,press_y2,press_y1];
            top_press_rect = [3,4,press_x3,press_x3,press_x4,press_x4,press_y3,press_y4,press_y4,press_y3];
            pad1 = zeros(1,4);
            pad2 = zeros(2,4);
            pad4 = zeros(4,4);
            % press_support poly
            % x11 = x7+self.support_poly_gap/sin(deg2rad(self.theta));
            sup_x2 = x7+self.support_poly_gap/sin(deg2rad(self.theta));
            sup_x1 = sup_x2+(self.angle_wall_len-self.support_poly_gap/sin(deg2rad(self.theta)))*cos(deg2rad(self.theta));
            sup_x3 = x10-self.support_poly_gap/sin(deg2rad(self.theta));
            sup_x4 = sup_x3-3e-3;
            sup_y1 = y1+poly_y_off+self.support_poly_gap;
            sup_y2 = y2;
            press_support_poly = [2,4,sup_x1,sup_x2,sup_x3,sup_x4,sup_y1,sup_y2,sup_y2,sup_y1];
            top_key_rect = [3,4,sup_x2,sup_x2,sup_x3,sup_x3,top_key_y1,top_key_y2,top_key_y2,top_key_y1];
            geo_matrix = [main_rect pad1;
                          corner_rects pad4;
                          left_poly;
                          right_poly;
                          mid_polys pad2;
                          press_support_poly pad1;
                          top_key_rect pad1;
                          bottom_press_rect pad1;
                          top_press_rect pad1]';
            [dl,bt] = decsg(geo_matrix, self.geo_formula, self.geo_names);
            [dl2,bg2] = csgdel(dl,bt);
            self.edges = dl2;
        end

        function setModel(self, force_app)
            self.force_app = force_app;
            self.surface_traction = self.force_app/self.total_thk/self.press_len;
            self.model = createpde(2);
            geometryFromEdges(self.model,self.edges);
            specifyCoefficients(self.model,'m',0,'d',0,'c',self.pde_c,'a',0,'f',self.pde_f);
            applyBoundaryCondition(self.model,'neumann','edge',1:self.model.Geometry.NumEdges,'g',[0 0]);
            fixed_edge = intersect(find(self.edges(4,:)==min(self.edges(4,:))),find(self.edges(5,:)==min(self.edges(5,:))));
            force_edge = intersect(find(self.edges(4,:)==max(self.edges(4,:))),find(self.edges(5,:)==max(self.edges(5,:))));
            applyBoundaryCondition(self.model,'dirichlet','edge',fixed_edge,'u',[0 0]);
            applyBoundaryCondition(self.model,'neumann','edge',force_edge,'g',[0 -self.surface_traction]);
            % setInitialConditions(self.model,0);
            generateMesh(self.model,'Hmax',.5*1e-3,'GeometricOrder','quadratic');
        end

        function solveModel(self)
            self.results = solvepde(self.model);
            [self.results_p,self.results_e,self.results_t]=meshToPet(self.results.Mesh);
            self.deflection_x = max(self.results.NodalSolution(:,1));
            self.deflection_y = abs(min(self.results.NodalSolution(:,2)));
        end

        function showGeometry(self)
            pdegplot(self.edges,'EdgeLabels','on','FaceLabels','on');
        end

        function showDeflectionX(self)
            pdeplot(self.results_p+self.results.NodalSolution',self.results_e,self.results_t,'XYData',self.results.NodalSolution(:,1));
        end

        function showDeflectionY(self)
            pdeplot(self.results_p+self.results.NodalSolution',self.results_e,self.results_t,'XYData',self.results.NodalSolution(:,2));
        end

        function calcDeflectionByForce(self, force_vector)
            self.force_app_buffer = force_vector;
            self.deflection_x_buffer = [];
            self.deflection_y_buffer = [];
            wall_thk = self.wall_thk*1e3;
            wall_len = self.wall_len*1e3;
            pinch_len = self.pinch_len*1e3;
            press_pos = self.press_pos*1e3;
            theta = self.theta;
            for force_app = self.force_app_buffer
                self.force_app = force_app;
                self.setAndSolve(wall_thk, wall_len, theta, force_app, pinch_len, press_pos);
                self.deflection_x_buffer = [self.deflection_x_buffer self.deflection_x];
                self.deflection_y_buffer = [self.deflection_y_buffer self.deflection_y];
            end
            self.deflection_x_buffer = reshape(self.deflection_x_buffer,[],length(self.force_app_buffer));
            self.deflection_y_buffer = reshape(self.deflection_y_buffer,[],length(self.force_app_buffer));
        end

        function calcDeflectionByForceAndAngle(self, angle_vector, force_app_vector)
            self.theta_buffer = angle_vector;
            self.force_app_buffer = force_app_vector;
            self.deflection_x_buffer = [];
            self.deflection_y_buffer = [];
            wall_thk = self.wall_thk*1e3;
            wall_len = self.wall_len*1e3;
            pinch_len = self.pinch_len*1e3;
            press_pos = self.press_pos*1e3;
            for force_app = self.force_app_buffer
                self.force_app = force_app;
                for theta = self.theta_buffer
                    self.setAndSolve(wall_thk, wall_len, theta, force_app, pinch_len, press_pos);
                    self.deflection_x_buffer = [self.deflection_x_buffer self.deflection_x];
                    self.deflection_y_buffer = [self.deflection_y_buffer self.deflection_y];
                end
            end
            self.deflection_x_buffer = reshape(self.deflection_x_buffer,[],length(self.force_app_buffer));
            self.deflection_y_buffer = reshape(self.deflection_y_buffer,[],length(self.force_app_buffer));
        end

        function calcDeflectionByForcePosition(self, pos_vector)
            self.force_pos_buffer = pos_vector;
            self.deflection_x_buffer = [];
            self.deflection_y_buffer = [];
            wall_thk = self.wall_thk*1e3;
            wall_len = self.wall_len*1e3;
            pinch_len = self.pinch_len*1e3;
            for press_pos = self.force_pos_buffer
                self.setAndSolve(wall_thk, wall_len, self.theta, self.force_app, pinch_len, press_pos);
                self.deflection_x_buffer = [self.deflection_x_buffer self.deflection_x];
                self.deflection_y_buffer = [self.deflection_y_buffer self.deflection_y];
            end
            self.deflection_x_buffer = reshape(self.deflection_x_buffer,[],length(self.force_pos_buffer));
            self.deflection_y_buffer = reshape(self.deflection_y_buffer,[],length(self.force_pos_buffer));
        end
    end
end