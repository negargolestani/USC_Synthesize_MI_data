classdef TWOPORTNET < handle
    % This class represents an two port network     
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        elements 
        f
        Z0
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function self = TWOPORTNET(elements, f, varargin)
            % Class constructor
            %       TWOPORTNET(elements, f)
            %       TWOPORTNET(elements, f, Z0)
            %
            % INPUTS:
            %       elements (N_by_1 cell array of elements): elements can
            %       be INDUCTOR, RESISTOR, CAPACITOR
            %       f (1_by_1 duoble): frequency (Hz)
            %       Z0 (1_by_1 double): charachteristic impedance (ohm)
            %                       
            % OUTPUTS:
            %       TWOPORTNET object
            
            self.elements = elements;                                      
            self.f = f;   
            if nargin == 3
                Z0 = varargin{1};
            else
                Z0 = 50;
            end
            self.Z0 = Z0;
        end
        % -----------------------------------------------------------------
        function [ A, B, C, D ] = ABCD(self)
            %  This function returns ABCD-parameters 
            %
            % INPUTS:
            %        None
            %
            % OUTPUTS:
            %       A, B, C, D  (1_by_1 double): ABCD parameters 
            
            ABCDmtx = eye(2);
            for n = 1:length(self.elements)
                ABCDmtx = ABCDmtx * self.elements{n}.ABCD(self.f);
            end
            A = ABCDmtx(1,1); B = ABCDmtx(1,2); 
            C = ABCDmtx(2,1); D = ABCDmtx(2,2);
        end
        % -----------------------------------------------------------------
        function [ Z11, Z21, Z12, Z22 ] = Z(self)
            %  This function returns Z-parameters  
            %
            % INPUTS:
            %        None
            %
            % OUTPUTS:
            %       Z11, Z21, Z12, Z22 (1_by_1 double): Z parameters 
            
            [A, B, C, D] = self.ABCD();
            Z11 = A ./ C;
            Z12 = ( A .* D - B .* C ) ./ C;
            Z21 = 1 ./ C;
            Z22 = D ./ C;            
        end
        % -----------------------------------------------------------------
        function [ S11, S21, S12, S22 ] = S(self)
            %  This function returns S-parameters  
            %
            % INPUTS:
            %        None
            %
            % OUTPUTS:
            %       S11, S21, S12, S22 (1_by_1 double): S parameters 
            
            [ A, B, C, D ] = self.ABCD();            
            S11 = ( A + B/self.Z0 - C*self.Z0 - D) ./ ( A + B/self.Z0 + C*self.Z0 + D);
            S12 = 2*( A.*D - B.*C ) ./ ( A + B/self.Z0 + C*self.Z0 + D);
            S21 = 2 ./ ( A + B/self.Z0 + C*self.Z0 + D);
            S22 = ( -A + B/self.Z0 - C*self.Z0 + D) ./ ( A + B/self.Z0 + C*self.Z0 + D);            
        end
        % -----------------------------------------------------------------
        function [ S11_dB, S21_dB, S12_dB, S22_dB ] = SdB(self)
            %  This function returns SdB-parameters 
            %
            % INPUTS:
            %        None
            %
            % OUTPUTS:
            %       S11_dB, S21_dB, S12_dB, S22_dB (1_by_1 double): 
            %       S(dB) parameters 
            
            [S11, S21, S12, S22] = self.S( self.Z0 );
            S11_dB = 20*log10( abs(S11)) + 1i* (angle(S11));
            S21_dB = 20*log10( abs(S21)) + 1i* (angle(S21));
            S12_dB = 20*log10( abs(S12)) + 1i* (angle(S12));
            S22_dB = 20*log10( abs(S22)) + 1i* (angle(S22));
        end
        % -----------------------------------------------------------------
    end
    
end

