classdef FEMComsolSolverClass
    %FEMCOMSOLSOLVERCLASS solutionIndex
    % is used to perfrom FEM simulation on a given geometry   
    % by creating,configuring COMSOL models and running simulations on a client server with the COMSOL Matlab livelink functionality.   
    
    properties
        model                   % COMSOL model object 
        modelName               % name of the COMSOL model object 
        geometryName            % name of the imported geometry 
        geometry                % Comsol geometry tag 
        componentName           % name of the COMSOL component within the model 
        volumes                 % geometryclass objects
    end
    %%%Additional properties for Comsol imported in defVolumes
    %   volumes.wp_geom           %xy/yz/zx WorkPlane geometry, 3rd coordinate inside detector_volume
    %   volumes.offsets           %WorkPlane 3rd coordinate offset
    %   volumes.sel               %sel index for Explicit Selection in Comsol
    %   volumes.selNames          %sel names from user
    %   volumes.mergeOperation    %Operation "-" is Boolean is used
    %   volumes.dif               %dif index for Comsol difference operation
    %   volumes.primitive         %Primitive name in Comsol used as an input for difference operation

    methods
        %% Constructor method 
        function obj = FEMComsolSolverClass()

        end
        %% Create a COMSOL model
        function obj = createModel(obj,modelName,componentName,geometryName)
            import com.comsol.model.util.*;
            ModelUtil.create(modelName); 
            obj.model = ModelUtil.model(modelName); 
            obj.modelName = modelName; 
            obj.model.component.create(componentName,false);
            obj.componentName = componentName;
            obj.geometryName = geometryName; 
            obj.geometry=obj.model.component(componentName).geom.create(geometryName,3);
                            
        end 

        function obj = importParams(obj, geometry)
            [~, selsize] = size(geometry.selections);
            for ss = 1:selsize
                [~, mergesize] = size(geometry.selections(ss).mergelist);
                for ms = 1:mergesize
                    [~, volumesize] = size(geometry.volumes);
                    for vs = 1:volumesize
                        if contains(geometry.selections(ss).mergelist{ms},geometry.volumes(vs).name)%contains pixel
                            ox = geometry.volumes(vs).origin(1);
                            obj.model.param.set(strcat(geometry.volumes(vs).name,'_ox'),ox);
                            oy = geometry.volumes(vs).origin(2);
                            obj.model.param.set(strcat(geometry.volumes(vs).name,'_oy'),oy);
                            oz = geometry.volumes(vs).origin(3);
                            obj.model.param.set(strcat(geometry.volumes(vs).name,'_oz'),oz);
                            vx = geometry.volumes(vs).vertexes(1);
                            obj.model.param.set(strcat(geometry.volumes(vs).name,'_vx'),vx);
                            vy = geometry.volumes(vs).vertexes(2);
                            obj.model.param.set(strcat(geometry.volumes(vs).name,'_vy'),vy);
                            vz = geometry.volumes(vs).vertexes(3);
                            obj.model.param.set(strcat(geometry.volumes(vs).name,'_vz'),vz);
                            if ~isempty(find([vx vy vz]==0))
                                idp = find([vx vy vz]==0);
                                Tbl = [ox oy oz];
                                offset = Tbl(idp);
                                obj = defVolumes(obj, geometry.volumes(vs), idp, offset,vs,geometry.selections(ss).name,geometry.selections(ss).mergelist{ms});
                            else
                                obj = defVolumes(obj, geometry.volumes(vs), NaN, NaN,vs,geometry.selections(ss).name,geometry.selections(ss).mergelist{ms});
                            end
                        end
                    end
                end
            end
        end

        function obj = defVolumes(obj,volume,idPlane,offset,volumeNum,selectionName,mergeOperation)
            volume.wp_geom=num2str(NaN);
            volume.offsets=num2str(offset);
            volume.sel=strcat('sel',num2str(volumeNum));
            volume.selName=selectionName;
            volume.mergeOperation=mergeOperation;
            volume.dif=[];
            volume.primitive=[];
            obj.volumes = cat(2,obj.volumes,volume);
            if idPlane==1
                obj.volumes(end).wp_geom='yz';
            elseif idPlane==2
                obj.volumes(end).wp_geom='zx';
            elseif idPlane==3
                obj.volumes(end).wp_geom='xy';
            end
        end
        
        %% Load a COMSOL model saved in local path   
        function obj = loadModel(obj,filePath,fileName)
             presentFolder = cd(filePath); 
             obj.model = mphopen(fileName); 
             cd(presentFolder)
        end 
        
        %% Save COMSOL model to a local path 
        function saveModel(obj,filePath,fileName)
            presentFolder = cd(filePath);
            mphsave(obj.model,fileName); 
            cd(presentFolder); 
        end 
        %% Import Geometry into COMSOL model 
        function obj = importGeometry(obj)
            g = obj.geometry;
            Num=size(obj.volumes,2);
            for l=1:Num
                if contains(obj.volumes(l).selName,"detector_volume")
                    component=obj.volumes(l).name;
                    g.create(component, 'Block');
                    g.feature(component).set('selresult', true);
                    g.feature(component).set('selresultshow', 'all');
                    g.feature(component).set('pos', {'detector_ox' 'detector_oy' 'detector_oz'});
                    g.feature(component).set('size', {'detector_vx' 'detector_vy' 'detector_vz'});

                    %%%Find workPlanes parameters
                    N=size(obj.volumes,2);
                    wpList={obj.volumes.wp_geom};
                    offsetList={obj.volumes.offsets};
                    pairs=string(zeros(1, N));
                    for i=1:N
                        pairs(i)=[wpList{i},offsetList{i}];
                        if contains(component,obj.volumes(i).name)
                            %sel=obj.volumes(i).sel;
                            obj.volumes(i).primitive=component;
                        end
                    end

                    %%% Remove NaN indexing
                    uniquePairs = unique(pairs,'stable');
                    idx_NaN = find(strcmp(uniquePairs, 'NaNNaN'));
                    uniquePairs(idx_NaN)=[];

                    %%%Create all unique pairs offset + WorkPlane geometry
                    for i=1:size(uniquePairs,2)
                        wp_geom=extract(uniquePairs(i),lettersPattern);
                        offset=str2double(erase(uniquePairs(i),wp_geom));
                        g.create(strcat('wp',num2str(i)), 'WorkPlane');
                        obj.model.component(obj.componentName).geom(obj.geometryName).feature(strcat('wp',num2str(i))).set('quickplane', wp_geom);
                        g.feature(strcat('wp',num2str(i))).set('quickz', offset);
                        g.feature(strcat('wp',num2str(i))).set('selplaneshow', true);
                        g.feature(strcat('wp',num2str(i))).set('unite', true);
                        g.feature(strcat('wp',num2str(i))).geom.create('cro1', 'CrossSection');
                        for j=1:size(obj.volumes,2)
                            if strcmp(obj.volumes(j).wp_geom,wp_geom) && str2double(obj.volumes(j).offsets)==offset
                                obj.volumes(j).wp=strcat('wp',num2str(i));
                            end
                        end
                    end
                end
            end

            for l=1:Num
                if contains(obj.volumes(l).selName,'electrode')
                    component=obj.volumes(l).name;
                    pix_id=str2double(extract([obj.volumes(l).sel],digitsPattern));
                    r_num=strcat('r',num2str(pix_id));
                    obj.volumes(l).primitive=r_num;
                    wp=obj.volumes(l).wp;
                    sel=obj.volumes(l).sel;
                    wp_geom=obj.volumes(l).wp_geom;

                    %%%Create rectangle for the pixel
                    g.feature(wp).geom.create(r_num, 'Rectangle');
                    g.feature(wp).geom.feature(r_num).label(component);
                    g.feature(wp).geom.feature(r_num).set('base', 'corner');

                    %%%Define the rectangle sizes
                    if strcmp(wp_geom,'yz')
                        g.feature(wp).geom.feature(r_num).set('pos', [strcat(component,'_oy') strcat(component,'_oz')]);
                        g.feature(wp).geom.feature(r_num).set('size', [strcat(component,'_vy') strcat(component,'_vz')]);
                    elseif strcmp(wp_geom,'zx')
                        g.feature(wp).geom.feature(r_num).set('pos', [strcat(component,'_oz') strcat(component,'_ox')]);
                        g.feature(wp).geom.feature(r_num).set('size', [strcat(component,'_vz') strcat(component,'_vx')]);
                    elseif strcmp(wp_geom,'xy')
                        g.feature(wp).geom.feature(r_num).set('pos', [strcat(component,'_ox') strcat(component,'_oy')]);
                        g.feature(wp).geom.feature(r_num).set('size', [strcat(component,'_vx') strcat(component,'_vy')]);
                    end

                    %%%Create explicit selection for the Pixel
                    g.feature(wp).geom.create(sel, 'ExplicitSelection');
                    g.feature(wp).geom.feature(sel).label(strcat('Pix_',num2str(pix_id)));
                    g.feature(wp).geom.feature(sel).selection('selection').set(strcat('r',num2str(pix_id),'(1)'), 1);
                end
            end

            for l=1:Num
                %%%Find mergelist - all the pixels in electrode
                mergelist=find(strcmp([obj.volumes(:).selName],obj.volumes(l).selName));
                %%%Find if the mergelist > 1 and there is a boolean operation
                if length(mergelist)>1 && any(contains([obj.volumes(mergelist).mergeOperation],'-'))
                    %%%Check if the boolean operation is in one WorkPlane
                    if all(strcmp({obj.volumes(mergelist).wp}, obj.volumes(mergelist(1)).wp))
                        %%%If the boolean operation is already made (One pixel can be assigned to one selection)
                        if isempty([obj.volumes(mergelist).dif])
                            %%%Find any difference operation maded previously, assign the next index
                            if any([obj.volumes(:).dif])
                                %%%Max number in differences definition
                                dif_id=max(str2double(extract([obj.volumes(:).dif],digitsPattern)))+1;
                            else
                                dif_id=1;
                            end
                            %%%Find which primitives to add and which primitives to subtract
                            idlogic_subtract=contains([obj.volumes(mergelist).mergeOperation],'-');
                            dif_num=strcat('dif',num2str(dif_id));

                            %%%Assign dif number for all the pixels within current dif operation
                            [obj.volumes(mergelist).dif]=deal(dif_num);
                            g.feature(obj.volumes(mergelist(1)).wp).geom.create(dif_num, 'Difference');

                            %%%Add primitives
                            g.feature(obj.volumes(mergelist(1)).wp).geom.feature(dif_num).selection('input').set({obj.volumes(mergelist(~idlogic_subtract)).primitive});
                            
                            %%%Subtract primitives
                            g.feature(obj.volumes(mergelist(1)).wp).geom.feature(dif_num).selection('input2').set({obj.volumes(mergelist(idlogic_subtract)).primitive});

                        end
                    end
                end
            end

            for l=1:Num
                if contains(obj.volumes(l).selName,"world_volume")
                    component=obj.volumes(l).name;
                    g.create(component, 'Block');
                    g.feature(component).set('selresult', true);
                    g.feature(component).set('selresultshow', 'all');
                    g.feature(component).set('base', 'corner');
                    g.feature(component).set('pos', [strcat(component,'_ox') strcat(component,'_oy') strcat(component,'_oz')]);
                    g.feature(component).set('size', [strcat(component,'_vx') strcat(component,'_vy') strcat(component,'_vz')]);

                    %%%Find 'world' selection
                    N=size(obj.volumes,2);
                    for i=1:N
                        if contains(component,obj.volumes(i).name)
                            [prim_wrld,obj.volumes(i).primitive]=deal(component);
                        elseif contains(obj.volumes(i).selName,"detector_volume")
                            prim_det=obj.volumes(i).primitive;
                        end
                    end

                    %%%Difference of detector and world. Final geometry
                    g.create('dif1', 'Difference');
                    g.feature('dif1').set('keepsubtract', true);
                    g.feature('dif1').selection('input').set(prim_wrld);
                    g.feature('dif1').selection('input2').set(prim_det);
                    g.run('fin');
                end
            end
        end

           
        %% Add materials setup simulation volume and mesh (simulation-ready)
        function obj = ComsolFEMInitialization(obj,detectorMaterial,epsr_det,worldMaterial,epsr_amb)

            %% Create materials for detector volume and world volume
            obj.model.component(obj.componentName).material.create(detectorMaterial, 'Common');
            obj.model.param.set('epsr_det',epsr_det);
            obj.model.component(obj.componentName).material.create(worldMaterial, 'Common');
            obj.model.param.set('epsr_amb',epsr_amb);

            %%%Assign relative permtitivity values
            obj.model.component(obj.componentName).material(detectorMaterial).propertyGroup('def').set('relativepermittivity', 'epsr_det');
            obj.model.component(obj.componentName).material(detectorMaterial).propertyGroup('def').set('relpermittivity', {'epsr_det' '0' '0' '0' 'epsr_det' '0' '0' '0' 'epsr_det'});
            obj.model.component(obj.componentName).material(worldMaterial).propertyGroup('def').set('relativepermittivity', 'epsr_amb');
            obj.model.component(obj.componentName).material(worldMaterial).propertyGroup('def').set('relpermittivity', {'epsr_amb' '0' '0' '0' 'epsr_amb' '0' '0' '0' 'epsr_amb'});

            %%%Assign material properties for detector_volume domain and world_volume domains
            obj.model.component(obj.componentName).material(worldMaterial).selection.set(1);
            obj.model.component(obj.componentName).material(detectorMaterial).selection.set(2);

            %% Declare physics interface
            obj.model.component(obj.componentName).physics.create('es1','Electrostatics',obj.geometryName);

            %%%All the terminal types - Voltage. Default voltage value is zero
            N=size(obj.volumes,2);
            for i=1:N
                if ~isempty(obj.volumes(i).wp)
                    term=strcat('term',num2str(i));
                    obj.model.component(obj.componentName).physics('es1').create(term, 'Terminal', 2);
                    obj.model.component(obj.componentName).physics('es1').feature(term).selection.named(strcat(obj.geometryName,'_',obj.volumes(i).wp,'_',obj.volumes(i).sel));
                    obj.model.component(obj.componentName).physics('es1').feature(term).set('TerminalType', 'Voltage');
                    obj.model.component(obj.componentName).physics('es1').feature(term).set('V0', 0);
                    obj.volumes(i).term=term;
                end
            end

           %% Generate mesh 
           obj.model.component(obj.componentName).mesh.create('mesh1','geom1'); 
           obj.model.component(obj.componentName).mesh.run
           
        end 
       
        %% Compute weighting field/potential for a given electrode 
        function obj = ComsolComputeWeigthingField(obj,electrodeName)

            %%% Declare and run the solution
            std_name = strcat('std1');
            std = obj.model.study.create(std_name);
            std.feature.create('stat', 'Stationary');
   
            N=size(obj.volumes,2);
            for i=1:N
                if strcmp(obj.volumes(i).selName,electrodeName)
                    obj.model.component(obj.componentName).physics('es1').feature(obj.volumes(i).term).set('V0', 1);
                end
            end

            std.run
        end

        %% Compute electric field given the voltages applied to the electrodes 
        function obj = ComsolComputeElectricField(obj,electrodeName,electrodeVoltage)

            %%% Declare and run the solution
            std_name = strcat('std1');
            std = obj.model.study.create(std_name);
            std.feature.create('stat', 'Stationary');

            N=size(obj.volumes,2);
            M=size(electrodeName);
            for j=1:M
                electrode=electrodeName(j);
                for i=1:N
                    if strcmp(obj.volumes(i).selName,electrode)
                        obj.model.component(obj.componentName).physics('es1').feature(obj.volumes(i).term).set('V0', electrodeVoltage(j));
                    end
                end
            end

            std.run
        end

    end

         
        
end

