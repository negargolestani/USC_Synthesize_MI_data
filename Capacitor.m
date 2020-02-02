classdef Capacitor
    % This class represents a capacitor
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        topology                                                            % toplogy of element
        value                                                               % value of element
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -----------------------------------------------------------------
        function self = Capacitor(topology, value)
            % Class constructor
            %       Capacitor(topology, value)
            %
            % INPUTS:
            %       topology (1_by_1 char): topology of element in circuit
            %                               'S': series
            %                               'P': parallel
            %       value (double): value of element (F)
            %
            % OUTPUTS:
            %       Capacitor object 
            
            self.topology = topology;
            self.value = value;
        end
        % -----------------------------------------------------------------
        function [ abcd ] = ABCD(self, f)
            % This function returns ABCD matrix of the element
            %
            % INPUTS:
            %       f (1_by_1 double): frequency (Hz)
            %
            % OUTPUTS:
            %       abcd (2_by_2 double): ABCD-parameters               
            
            w = 2*pi*f;
            z = 1./( 1i * w .* self.value);
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

