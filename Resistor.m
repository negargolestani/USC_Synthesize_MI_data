classdef Resistor
    % This class represents a resistor
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        topology                                                            % Toplogy of element
        value                                                               % value of element
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -----------------------------------------------------------------
        function self = Resistor(topology, value)
            % Class constructor
            %       Resistor(topology, value)
            %
            % INPUTS:
            %       topology (1_by_1 char): topology of element in circuit
            %                               'S': series
            %                               'P': parallel
            %       value (1_by_1 double): value of element (ohm)
            %
            % OUTPUTS:
            %       Resistor object
            
            self.topology = topology;
            self.value = value;
        end
        % -----------------------------------------------------------------
        function [ abcd ] = ABCD(self, varargin)
            % This function returns ABCD matrix of the element
            %
            % INPUTS:
            %       f (1_by_1 double): frequency (Hz)
            %
            % OUTPUTS:
            %       abcd (2_by_2 double): ABCD-parameters 
            
            z = self.value;
            switch self.topology
                case 'S'
                    abcd = [1, z; 0, 1];
                case 'P'
                    abcd = [1, 0; 1/z, 1];
            end
        end
        % -----------------------------------------------------------------
    end
end

