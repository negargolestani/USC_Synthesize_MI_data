classdef TRANSCEIVER < handle
    % This class defines MI transceiver with coil and matching circuit
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        coil                                                                % coil of transceiver
        matching = []                                                       % elements of matching circuit
        abcd
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function self = TRANSCEIVER(coil, varargin)
            % Class constructor
            %        TRANSCEIVER(coil, matchingE_1, matchingE_2,...)
            %
            % INPUTS:
            %       coil (1_by_1 Coil obj): coil of MI transceiver
            %       matchingE_i (1_by_1 <circuit element>):circuit elements
            %       of matching network
            %
            % OUTPUTS:
            %       TRANSCEIVER object
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
            %       newself (1_by_1 TRANSCEIVER obj): copy of coil object
            %       with all properties
            
            newself = TRANSCEIVER(self.coil);
            newself.matching = self.matching;
        end
        % -----------------------------------------------------------------
        function [ ] = move(self, new_loc, new_align)
            % This function moves/rotates TRANSCEIVER's coil to the given
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
        function [ twoPortModel ] = model(self, f)
            % This function returns two port circuit model of TRX
            % NOTE: No mutual inductance
            %
            % INPUTS:
            %       f (1_by_f double): frequency (Hz)
            %
            % OUTPUTS:
            %       elements (N_by_1 cell array): all elements of TRX
            
            coil_modelElements = self.coil.getmodelelements(f);
            if isempty(self.matching)
                elements = { coil_modelElements{end:-1:1} };
            else
                elements = { ...
                    coil_modelElements{end:-1:1}, ...
                    self.matching{end:-1:1} };
            end
            twoPortModel = TWOPORTNET(elements, f);
            [A, B, C, D] = twoPortModel.ABCD();                                       
            self.abcd = [A, B; C, D];                                       % also save ABCD parameters to seepd up calculations            
        end
        % -----------------------------------------------------------------
        
    end
    
end

