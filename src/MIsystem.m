classdef MISYSTEM < handle
    % This class represents an MI system with one reciever (RX) and
    % multiple transmitters (TXs)
    
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        f                                                                   % operating frequency
        Z0                                                                  % charachteristic Impedance
        RX                                                                  % reciever
        TX = []                                                             % transmitter
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % -----------------------------------------------------------------
        function self = MISYSTEM(f, RX, varargin)
            % Class constructor
            %       MISYSTEM(f, RX)
            %       MISYSTEM(f, RX, TX_1, TX_2, ...)
            %       MISYSTEM(f, RX, TX, Ntx)
            %
            % INPUTS:
            %       f (1_by_1 double): frequency (Hz)
            %       RX (1_by_1 Transceiver obj): receiver
            %       TX_n (1_by_1 Transceiver obj): transmitter
            %       Ntx (1_by_1 double): number of transmitters
            %
            % OUTPUTS:
            %       MISYSTEM object
            
            self.Z0 = 50;                                                    % charachteristic impedance
            self.f = f;                                                      % operating frequency of MI system
            self.setrx(RX);                                                  % reciever
            self.settx(varargin{:});                                        % transmitters
        end
        % -----------------------------------------------------------------
        function setrx(self, RX)
            % This function sets RX (transceiver) of MI system
            %
            % INPUTS:
            %        RX (1_by_1 Transceiver obj): receiver
            %
            % OUTPUTS:
            %       None
            
            self.RX = copy(RX);
            self.RX.model(self.f);
        end
        % -----------------------------------------------------------------
        function settx(self, varargin)
            % This function sets TXs (transceiver) of MI system
            %
            % INPUTS:
            %        TX_i (1_by_1 Transceiver obj): transmitter
            %        Ntx: number of transmitters
            %
            % OUTPUTS:
            %       None
            
            if nargin==3 && isnumeric(varargin{2})                          % Set same coil as TXs:  settx(Transceiver, 8)
                TXs = repelem(varargin{1},varargin{2});
            else                                                            % Set different coils as TXs: settx(Transceiver_1, Transceiver_2, ...)
                TXs = [varargin{:}];
            end
            
            for n = 1:length(TXs)
                self.TX = [self.TX, copy(TXs(n)) ];
                self.TX(n).model(self.f);
            end
            
        end                    
        % -----------------------------------------------------------------
        function [ synthData ] = synthesizedata(self, motions)
            % This function synthesizes S21_dB of MI system with given
            % motion data
            %
            % INPUTS:
            %        motions (1_by_Ncoils cell array of timeseries 
            %        collection):contain Location of coils' center and 
            %        Normal of the coil's surface 
            %        NOTE: first coil (motions{1}) contains loc/aligns of RX
            %
            % OUTPUTS:
            %        synthMI (1_by_Ntx cell array of timeseries): contains synthesize
            %        MI data generated correspond to each TX coil (S21(dB))            
                        
            Ntx = length(self.TX);
            time = motions{1}.Time;            
            Ntime = motions{1}.Length;
            S21_dB = zeros(Ntime, Ntx);
            
            for t = 1:Ntime
                
                % Move RX
                Crx = getdatasamples(motions{1}.Location, t);
                nrx = getdatasamples(motions{1}.Normal, t);
                self.RX.move (Crx, nrx);                                    
                
                for n = 1:Ntx
                    
                     % Move TXs
                    Ctx = getdatasamples(motions{n+1}.Location, t);
                    ntx = getdatasamples(motions{n+1}.Normal, t);
                    self.TX(n).move (Ctx, ntx);                            
                    
                    % Synthetic Data
                    M = self.RX.coil.mutualInductance(...
                        self.TX(n).coil, self.f);
                    ms_abcd = INDUCTOR('S',-M).ABCD(self.f);
                    mp_abcd = INDUCTOR('P', M).ABCD(self.f);
                    
                    % ABCD parameters -> S21 -> S21_dB
                    abcd = ...
                        self.RX.abcd * ...
                        ms_abcd * mp_abcd * ms_abcd * ...
                        self.TX(n).abcd;                    
                    A = abcd(1,1); B = abcd(1,2); 
                    C = abcd(2,1); D = abcd(2,2);
                    S21 = 2 ./ ( A + B/self.Z0 + C*self.Z0 + D);
                    S21_dB(t,n) = 20*log10( abs(S21)) + 1i* (angle(S21));
                end
            end
            
            % collection of timeseries object
            synthData = cell(Ntx,1);
            for n = 1:Ntx
                synthData{n} = timeseries( ...
                    fillmissing( S21_dB(:,n), 'linear'), ...                % Clean data: remove NaNs
                    time, ...
                    'Name', ['TX',num2str(n)] );
                synthData{n}.DataInfo.Units = 'dB';
            end
            synthData = tscollection(synthData, 'Name','Synthetic MI Data');
        end
        % -----------------------------------------------------------------        
    end
end

