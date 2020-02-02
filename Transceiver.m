classdef Transceiver < handle
    % This class defines MI transceiver with coil and matching circuit    
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        coil                                                                % coil of transceiver 
        matching = []                                                       % elements of matching circuit
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function self = Transceiver(coil, varargin)
            % Class constructor
            %        Transceiver(coil, matchingE_1, matchingE_2,...)
            %
            % INPUTS:
            %       coil (1_by_1 Coil obj): coil of MI transceiver 
            %       matchingE_i (1_by_1 <circuit element>):circuit elements
            %       of matching network
            %
            % OUTPUTS:
            %       Transceiver object   
            %
            % NOTE:
            %       <circuit element>: Resistor, Capacitor, Inductor 
            
            if nargin >= 1
                self.coil = copy(coil);                
                if nargin > 1
                    self.matching = varargin;
                end
            end
        end
        % -----------------------------------------------------------------
        function [ newself ]  = copy(self)
            % This function makes a copy of object.
            %
            % INPUTS:
            %      None  
            %
            % OUTPUTS:
            %       newself (1_by_1 Transceiver obj): copy of coil object 
            %       with all properties
            
            newself = Transceiver(self.coil);
            newself.matching = self.matching;
        end
        % -----------------------------------------------------------------
        function [ ] = move(self, new_loc, new_align)
            % This function moves/rotates Transceiver's coil to the given 
            % location/alignment
            %
            % INPUTS:
            %       new_loc (1_by_3 vector): coils' center location (m)
            %       new_align (1_by_3 vector): coils' normal surface 
            %
            % OUTPUTS:
            %       None     
            
            self.coil.move(new_loc, new_align);            
        end
        % -----------------------------------------------------------------
        
    end
    
end

