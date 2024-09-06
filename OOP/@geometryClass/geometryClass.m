classdef geometryClass
    %
    %
    %DESCRIPTION: 
    %It defines the detector geometry has an ensemble of cuboids, grouped in
    %3 major functional groups:
    %   -detector volume
    %   -electrodes 
    %   -world volume 
    
    properties
        volumes         % Array of volumes representing the detector geoemtry  
        selections         % world volume, detector volume or selections 
    end
    
    methods
        % Constructor method
        function obj = geometryClass()       
            obj.volumes = []; 
            obj.selections = []; 
        end

        obj = buildPixelDetector(obj, pitch, arrayDims, interpixelGap, guardDepth, thickness, worldOffset)
        

   
        function obj = addPrimitive(varargin)
            
           obj = varargin{1}; 
           
           if varargin{2} == "primitive"
               volume.name = varargin{3}; 
               volume.type = varargin{2}; 
               volume.origin = varargin{4}; 
               volume.vertexes = varargin{5};
               obj.volumes = cat(2,obj.volumes,volume);  
           end 
           
        end
        
        function obj = definePlane(varargin)
            
           obj = varargin{1}; 
           
           if varargin{2} == "plane"
               volume.name = varargin{3}; 
               volume.type = varargin{2}; 
               volume.origin = varargin{4}; 
               volume.vertexes = varargin{5};
               obj.volumes = cat(2,obj.volumes,volume);  
           end 
           
        end
        % Merge volumes 
        function obj = addSelections(varargin)
            % Merges volumes into functional groups.
            %
            % Parameters:
            %   varargin: Variable input arguments.
            %       - varargin{1}: Instance of geometryClass.
            %       - varargin{2}: Name of the functional group.
            %       - varargin{3:end}: Volumes to be merged.
            %
            % Returns:
            %   obj: Updated instance of geometryClass.

            obj = varargin{1}; 
            [~, listsize] = size(obj.selections);   
            obj.selections(listsize+1).name = varargin{2}; 
            obj.selections(listsize+1).type = "selection";           
            obj.selections(listsize+1).mergelist = varargin(3:end);  
        end 
        

        function plotVolumes(obj)
            % Plots all volumes in the geometry.
            
            figure("Name", "Detector Geometry")
            
            for nv = 1 : length(obj.volumes)
                if obj.volumes(nv).type == "primitive"
                    plotcube(obj,obj.volumes(nv).vertexes,obj.volumes(nv).origin,.1,[1 0 0]); 
                end 
            end
            xlabel("x distance (m)")
            ylabel("y distance (m)")
            zlabel("z distance (m)")
        end
            
        %--------------------other functions-------------------%
        plotcube(obj,varargin)

       
    end
end

