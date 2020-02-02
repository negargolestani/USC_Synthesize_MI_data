classdef MIsystem < handle
    % This class represents an MI system with one reciever (RX) and 
    % multiple transmitters (TXs)

    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f                                                                   % operating frequency
        Z0                                                                  % charachteristic Impedance         
        RX                                                                  % reciever
        TX = []                                                             % transmitter
        
        % These parameters are used to speed up the calculation
        RX_abcd                                                             % ABCD-parameters of RX
        TX_abcd = []                                                        % ABCD-parameters of TXs
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function self = MIsystem(f, RX, varargin)
            % Class constructor
            %       Coil(f, RX)
            %       Coil(f, RX, TX_1, TX_2, ...)
            %       Coil(f, RX, TX, Ntx)
            %
            % INPUTS:
            %       f (1_by_1 double): frequency (Hz)
            %       RX (1_by_1 Transceiver obj): receiver
            %       TX_n (1_by_1 Transceiver obj): transmitter
            %       Ntx (1_by_1 double): number of transmitters 
            %
            % OUTPUTS:
            %       MIsystem object 
            
            self.Z0 = 50;                                                    % charachteristic impedance 
            self.f = f;                                                      % operating frequency of MI system
            self.set_RX(RX);                                                 % reciever
            self.set_TX(varargin{:});                                        % transmitters
            self.save_ABCD();                                                % save ABCD parameters of RX/TXs
        end
        % -----------------------------------------------------------------
        function set_RX(self, RX)
            % This function sets RX (transceiver) of MI system
            %
            % INPUTS:
            %        RX (1_by_1 Transceiver obj): receiver
            %
            % OUTPUTS:
            %       None
            
            self.RX = copy(RX);
        end
        % -----------------------------------------------------------------
        function set_TX(self, varargin)
            % This function sets TXs (transceiver) of MI system
            %
            % INPUTS:
            %        TX_i (1_by_1 Transceiver obj): transmitter
            %        Ntx: number of transmitters             
            %
            % OUTPUTS:
            %       None       
            
            if nargin==3 && isnumeric(varargin{2})                          % Set same coil as TXs:  set_TX(Transceiver, 8)
                TXs = repelem(varargin{1},varargin{2});
            else                                                            % Set different coils as TXs: set_TX(Transceiver_1, Transceiver_2, ...)
                TXs = [varargin{:}];
            end
            
            for n = 1:length(TXs)
                self.TX = [self.TX, copy(TXs(n)) ];
            end
            
        end
        % -----------------------------------------------------------------        
        function [ ] = save_ABCD(self)
            % This function saves two-port network model of RX and TXs to
            % speed up calculation
            %
            % INPUTS:
            %        None             
            %
            % OUTPUTS:
            %       None  
            
            RX_elements = [ ...
                self.RX.matching, ...
                self.RX.coil.get_modelElements(self.f)...
                ] ;
            self.RX_abcd = eye(2);
            for e = 1:length(RX_elements)
                self.RX_abcd = self.RX_abcd * RX_elements{e}.ABCD(self.f);
            end
            
            for n = 1:length(self.TX)
                coil_modelElements = self.TX(n).coil.get_modelElements(self.f);
                
                if isempty(self.TX(n).matching)
                    TX_elements = { coil_modelElements{end:-1:1} };
                else                    
                    TX_elements = { ...
                        coil_modelElements{end:-1:1}, ...
                        self.TX(n).matching{end:-1:1} };                   
                end    
                self.TX_abcd{n} = eye(2);
                for e = 1:length(TX_elements)
                    self.TX_abcd{n} = self.TX_abcd{n} * TX_elements{e}.ABCD(self.f);
                end
            end
        end
        % -----------------------------------------------------------------
        function [ A, B, C, D ] = ABCD(self)
            %  This function returns ABCD-parameters between RX and each TX 
            %
            % INPUTS:
            %        None             
            %
            % OUTPUTS:
            %       A, B, C, D (each one is 1_by_Ntx double): 
            %       ABCD parameters of MI system 
                    
            for n = 1:length(self.TX)
                M = self.RX.coil.mutualInductance(self.TX(n).coil, self.f);
                ms_abcd = Inductor('S',-M).ABCD(self.f);
                mp_abcd = Inductor('P', M).ABCD(self.f);
                
                abcd = ...
                    self.RX_abcd * ...
                    ms_abcd * mp_abcd * ms_abcd * ...
                    self.TX_abcd{n};
                
                A(n) = abcd(1,1);
                B(n) = abcd(1,2);
                C(n) = abcd(2,1);
                D(n) = abcd(2,2);
            end
        end
        % -----------------------------------------------------------------
        function [ Z11, Z21, Z12, Z22 ] = Z(self)
            %  This function returns Z-parameters between RX and each TX 
            %
            % INPUTS:
            %        None             
            %
            % OUTPUTS:
            %       Z11, Z21, Z12, Z22 (each one is 1_by_Ntx double): 
            %       Z-parameters of MI system 
            %                       
            
            [A, B, C, D] = self.ABCD( );
            Z11 = A ./ C;
            Z12 = ( A .* D - B .* C ) ./ C;
            Z21 = 1 ./ C;
            Z22 = D ./ C;
        end
        % -----------------------------------------------------------------
        function [ S11, S21, S12, S22 ] = S(self)
            %  This function returns S-parameters between RX and each TX 
            %
            % INPUTS:
            %        None             
            %
            % OUTPUTS:
            %       A, B, C, D (each one is 1_by_Ntx double): 
            %       S-parameters of MI system 
            
            [A, B, C, D] = self.ABCD();
            S11 = ( A + B/self.Z0 - C*self.Z0 - D) ./ ( A + B/self.Z0 + C*self.Z0 + D);
            S12 = 2*( A.*D - B.*C ) ./ ( A + B/self.Z0 + C*self.Z0 + D);
            S21 = 2 ./ ( A + B/self.Z0 + C*self.Z0 + D);
            S22 = ( -A + B/self.Z0 - C*self.Z0 + D) ./ ( A + B/self.Z0 + C*self.Z0 + D);
        end
        % -----------------------------------------------------------------
        function [ S11_dB, S21_dB, S12_dB, S22_dB ] = SdB(self)
            %  This function returns SdB-parameters between RX and each TX 
            %
            % INPUTS:
            %        None             
            %
            % OUTPUTS:
            % OUTPUTS:
            %       A, B, C, D (each one is 1_by_Ntx double): 
            %       SdB-parameters of MI system )  
            
            [S11, S21, S12, S22] = self.S( );
            S11_dB = 20*log10( abs(S11)) + 1i* (angle(S11));
            S21_dB = 20*log10( abs(S21)) + 1i* (angle(S21));
            S12_dB = 20*log10( abs(S12)) + 1i* (angle(S12));
            S22_dB = 20*log10( abs(S22)) + 1i* (angle(S22));
        end
        % -----------------------------------------------------------------  
                
    end
    
end

