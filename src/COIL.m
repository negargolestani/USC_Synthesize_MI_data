classdef COIL < handle
    % This class represents an air-cored COILs with circular cross-section
    % wire
      
    properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        a                                                                   % coil radius 
        b                                                                   % coil length
        N                                                                   % number of turns 
        phi_w                                                               % coil wire diameter 
        loc                                                                 % location of coil's center 
        align                                                               % normlaized alignment of coils' surface normal 
    end
    
    methods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % -----------------------------------------------------------------
        function self  = COIL(a, N, awg, varargin)
            % Class constructor
            %               COIL(a, N, awg)
            %               COIL(a, N, awg, loc, align)
            %
            % INPUTS:
            %        a (1_by_1 double): coil radius (m)
            %        N (1_by_1 double): coil Number of turns
            %        awg (1_by_1 double): coil's wire AWG
            %
            %        loc (3_by_1 double): coil's center location (m)
            %        align (3_by_1 double): coil's surface normal alignment  
            %
            % OUTPUTS:
            %       COIL object
            
            if nargin >= 3
                PHI_W = 1e-3*[ ...
                    8.2515, 7.3481, 6.5437, 5.8273, 5.1894, 4.6213, ...
                    4.1154,3.6649, 3.2636, 2.9064, 2.5882, 2.3048, ...
                    2.0525, 1.8278, 1.6277, 1.4495, 1.2908, 1.1495, ...
                    1.0237, 0.9116, 0.8118, 0.7229, 0.6438, 0.5733, ...
                    0.5106, 0.4547, 0.4049, 0.3606, 0.3211, 0.2859, ...
                    0.2546, 0.2268, 0.2019, 0.1798, 0.1601, 0.1426, ...
                    0.1270, 0.1131, 0.1007, 0.0897,0.0799, 0.071, ...
                    0.064, 0.056, 0.051, 0.045,0.040, 0.036, 0.031, ...
                    0.028,0.025];
                self.a = a;
                self.N = N;
                self.phi_w = PHI_W(awg+1);
                self.b =  self.N*2*self.phi_w;
                
                if nargin == 5
                    loc = varargin{1};
                    align = varargin{2};
                    align = align ./ sqrt(sum(align.^2));                   % Normalized alignment
                    
                    self.loc = loc;
                    self.align = align ./ sqrt(sum(align.^2));
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
            %       newself (COIL obj): copy of coil with all properties
            
            newself = feval(class(self));                                   % instantiate new object of the same class.
            
            p = properties(self);                                           % copy all properties.
            for i = 1:length(p)
                newself.(p{i}) = self.(p{i});
            end
        end
        % ----------------------------------------------------------------
        function [ Reff ] = effResistance(self, f)
            % This function returns effective resistance of coil
            %
            % INPUTS:
            %        f (1_by_1 double): frequency (Hz)
            %
            % OUTPUTS:
            %        Reff (1_by_1 double): effective resistance of coil (ohm)
            
            rho_w = 1.724e-8;                                               % copper wire resistivity
            mu0 = 4*pi*1e-7;                                                % permeability of free space
            e0 = 8.85418782e-12;                                            % electric constant of free space
            c = 299792458;                                                  % speed of light
            Z0 = sqrt(mu0/e0);                                              % impedance of free space
            lambda = c./f;
            beta = 2*pi ./ lambda;
            
            pitch = self.b ./self.N;
            l_w =  self.N .* sqrt(pi^2*(2*self.a).^2+pitch.^2);             % length of wire
            % l_w =  N .* 2*pi*a;                                           % length of wire
            delta_w = sqrt( rho_w ./ (pi * mu0 * f) );                      % skin depth
            
            Rrad = pi*self.N.^2*Z0.*(beta.*self.a).^4/6;
            Rloss = l_w .* rho_w./(pi.*self.phi_w.*delta_w);
            Reff = Rrad + Rloss;
            
            % Reff = (self.phi_w < delta_w).* (4 * rho_w * l_w ./ ...
            %    ( pi * self.phi_w.^2 )) + (self.phi_w >= delta_w).* ...
            %    (rho_w * l_w ./(pi*delta_w .* (delta_w) ));                % Effective Resistance
            % Reff = l_w * rho_w ./ ( pi * (self.phi_w/2).^2) * ...
            %    ( 1 + 1/48 * (self.phi_w/2 ./ delta_w).^2);
            
        end
        % ----------------------------------------------------------------
        function [ L ] = selfInductance(self)
            % This function returns self-inductance of coil
            %
            % OUTPUTS:
            %       None
            % 
            % OUTPUTS:
            %       L (1_by_1 double): self inductance of coil (H)
            
            mu0 = 4*pi*1e-7;
            L = mu0*self.a*self.N.^2 .* ...
                ( log(8*self.a./self.b) - 1/2 ...       
                + (self.b./self.a).^2/32 .* ...
                ( log(8*self.a./self.b)+1/4 )...
                - (self.b./self.a).^4/1024 .* ...
                ( log(8*self.a./self.b)-2/3 ) ...
                + (self.b./self.a).^6.*10/131072 .* ...
                ( log(8*self.a./self.b)-109/120 ) ...
                - (self.b./self.a).^8*35/4194304 .* ...
                ( log(8*self.a./self.b)-431/420 ) );                        % Coffin's Formula (for b~a)
            
%             L = (self.a*self.N)^2 / (22.9*self.a + 25.4*self.b)*1e-4 ;      % from River Publishers Series in Communications : Principles of Inductive Near Field Communications for Internet of Things
%             L = mu0 * self.N^2 * self.a * ...
%                 ( (1+self.phi_w^2/8/self.a^2)* ...   
%                log(8*self.a/self.phi_w) ...
%                + self.phi_w^2/24/self.a^2 - 1.75 );                         % Rayleigh and Niven's formula
        end
        % ----------------------------------------------------------------
        function [ M ] = mutualInductance(self, coil_2, f)
            % This function returns mutual inductance between this coil and
            % given secondary coil
            %
            % INPUTS:
            %       coil_2 (COIL obj): second coil 
            %        f (1_by_1 double): frequency (Hz)
            %
            % OUTPUTS:
            %       M (1_by_1 double): mutual inductance between coil (H)
            
            [gamma, ~, ~] = dielectric( f , 'Air' );
            
            mu0 = 4*pi*1e-7;                                                % permeability of free space
            S_TX = pi * self.a^2;                                           % surface area of coils
            
            [coilsDistance, xRotAngle] = calculate_params( ...
                self.loc, self.align, coil_2.loc, coil_2.align );
            
            % Required parameters to calculate M
            a1 = coilsDistance(1);
            a2 = coilsDistance(2)*cos(xRotAngle) + coilsDistance(3)*sin(xRotAngle);
            B = sum(coilsDistance.^2);
            C =  coilsDistance(3)*cos(xRotAngle).*sin(xRotAngle) - ...
                coilsDistance(2)*sin(xRotAngle).^2 ;
            D =  coilsDistance(3)^2*cos(xRotAngle) - ...
                coilsDistance(2)*coilsDistance(3)*sin(xRotAngle);
            E = -cos(xRotAngle);
            f1 = -2*coilsDistance(1)*cos(xRotAngle);
            f2 = -coilsDistance(2)*cos(xRotAngle).^2 - coilsDistance(2) -...
                coilsDistance(3)*cos(xRotAngle).*sin(xRotAngle);
            G = -coilsDistance(1)^2*cos(xRotAngle)- coilsDistance(2)^2*cos(xRotAngle)-...
                coilsDistance(2)*coilsDistance(3)*sin(xRotAngle);
            
            fun =  @(rho,phi) ...
                1 ./ (2*pi*sqrt( rho.^2 + rho.*2.*(a1.*cos(phi)+a2.*sin(phi)) + B ).^5) .* ...
                real( exp(-gamma.*sqrt( rho.^2 + rho.*2.*(a1.*cos(phi)+a2.*sin(phi)) + B )) .* ...
                ( 1 + gamma.*sqrt( rho.^2 + rho.*2.*(a1.*cos(phi)+a2.*sin(phi)) + B ) ) ).* ...
                (C.*sin(phi).*rho + D) .* rho  + ...
                1 ./ (4*pi*sqrt( rho.^2 + rho.*2.*(a1.*cos(phi)+a2.*sin(phi)) + B ).^5) .* ...
                real( exp(-gamma.*sqrt( rho.^2 + rho.*2.*(a1.*cos(phi)+a2.*sin(phi)) + B )) .*...
                ( 1 + gamma.*sqrt( rho.^2 + rho.*2.*(a1.*cos(phi)+a2.*sin(phi)) + B )+ ...
                (gamma.*sqrt( rho.^2 + rho.*2.*(a1.*cos(phi)+a2.*sin(phi)) + B )).^2 ) ) .* ...
                (E*rho.^2 + (f1*cos(phi)+f2*sin(phi)).*rho + G ) .* rho ;
            M = mu0 * self.N * coil_2.N * S_TX * abs(integral2(fun,0,coil_2.a,-pi,pi));
            
        end
        % ----------------------------------------------------------------
        function [ ] = move( self, new_loc, new_align)
            % This function translates/rotates coil to the given 
            % location/alignment
            %
            % INPUTS:
            %       new_loc (1_by_3 vector): coils' center location (m)
            %       new_align (1_by_3 vector): coils' normal surface 
            %
            % OUTPUTS:
            %       None
            
            self.loc = new_loc;
            self.align = new_align ./ sqrt(sum(new_align.^2));
        end
        % ----------------------------------------------------------------
        function [ modelElements ] = getmodelelements(self, f)
            % This function returns elements of coils' circuit model
            %
            % INPUTS:
            %        f (1_by_1 double): frequency (Hz)
            %
            % OUTPUTS:
            %       modelElements (1_by_Nelements cell of <circuit element>)  
            %       all elements of circuit model
            %
            % NOTE:
            %       <circuit element>: RESISTOR, CAPACITOR, INDUCTOR             
            
            modelElements = {...
                RESISTOR('S', self.effResistance(f) ) ,...
                INDUCTOR('S', self.selfInductance() )...
                };
        end
        % ----------------------------------------------------------------
    end
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ coilsDistance, xRotAngle ] = calculate_params( loc_1, align_1, loc_2, align_2)
% This function calculates coil's distance and misalignment (xRotAngle)
% from coil's location and alignments for M calculation explained in paper* 
% loc_1 -> (0,0,0)  /	align_1 -> (0,0,1)  /	align_2 -> (0,ny,nz0) 
%
% INPUTS: 
%        loc_1 (1_by_3 double): first coil's center location (m)
%        align_1 (1_by_3 double): first coil's surface normal alignment
%        loc_2 (1_by_3 double): second coil's center location (m)
%        align_2 (1_by_3 double): second coil's surface normal alignment
%
% OUTPUTS:
%       coilsDistance (1_by_3 double): distance between coils (m)
%       xRotAngle (1_by_1 double): misalignmnet between coils (radians)


if align_1(2) == 0
    thetaX = 0;
else
    thetaX = atan( align_1(2)/align_1(3) );
end
thetaY = atan( -align_1(1)/ sqrt(align_1(2)^2+align_1(3)^2) );
align_2_new = get_rotationMatrix ( thetaX, thetaY , 0) * reshape(align_2,[3,1]);
if align_2_new(1) == 0
    thetaZ = 0;
else
    thetaZ = atan(align_2_new(1)/align_2_new(2));
end

Rtot = get_rotationMatrix ( thetaX, thetaY , thetaZ);
loc_2_new = Rtot * reshape ( loc_2-loc_1, [3,1]);
align_2_new = Rtot * reshape ( align_2, [3,1]);

coilsDistance = loc_2_new;
xRotAngle = atan(-align_2_new(2)/align_2_new(3));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ Rtotal ] = get_rotationMatrix(XrotAngle, YrotAngle, ZrotAngle)
% This function returns rotation matrix for Xrot, Yrot, Zrot
%
% INPUTS:
%       XrotAngle (1_by_1 double): rotation around x-axis (radians)
%       YrotAngle (1_by_1 double): rotation around y-axis (radians)
%       ZrotAngle (1_by_1 double): rotation around z-axis (radians)
% 
% OUTPUTS:
%       Rtotal (3_by_3 double): rotation around x,y,z axis in total 
%       Rtotal = Rz(theta_z) *  Ry(theta_y) *  Rx(theta_x) 


Rx = [ 1                 0                 0                 ; ...
       0                 cos(XrotAngle)    -sin(XrotAngle)   ; ...
       0                 sin(XrotAngle)    cos(XrotAngle)]   ;
   
Ry = [cos(YrotAngle)     0                 sin(YrotAngle)    ; ...
       0                 1                 0                 ; ...
      -sin(YrotAngle)    0                 cos(YrotAngle)]   ;
   
Rz = [cos(ZrotAngle)     -sin(ZrotAngle)   0                 ; ...
      sin(ZrotAngle)     cos(ZrotAngle)    0                 ; ...
      0                  0                 1]                ;


Rtotal = Rz*Ry*Rx;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ gamma, sigma, eps ] = dielectric(f , medium)
% This function returns dielectric constants of given medium 
% 
% INPUTS 
%       f (1_by_N double): freuency range (Hz)
%       medium (1_by_N char): 'Air' or human body tissue name 
% 
% OUTPUTS
%   gamma (1_by_N vector): complex propagation constant of medium
%   sigma (1_by_N vector): total conductivity of medium (S/m)
%   eps (1_by_N vector): real part of mediums' complex permittivity (F/m)

e0=8.854187817e-12;
mu0 = 4*pi*1e-7;

if strcmp(medium,'Air')
    gamma = 1i*2*pi*f*sqrt(mu0*e0);
    sigma = ones(size(f))*0;
    eps = ones(size(f))* e0;
else
    [gamma, sigma, eps] = dielectric_HB ( f, medium);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ gamma, sigma, eps ] = dielectric_HB(f, tissue)
% This function Returns Dielectric constants of Human body tissues
% Reference:
%	"Gabriel, Sami, R. W. Lau, and Camelia Gabriel. "The dielectric 
%   properties of biological tissues: II. Measurements in the frequency
%   range 10 Hz to 20 GHz." Physics in medicine & biology 41.11
%   (1996): 2251."
%
%   "Gabriel, Sami, R. W. Lau, and Camelia Gabriel. "The dielectric 
%   properties of biological tissues: III. Parametric models for the 
%   dielectric spectrum of tissues." Physics in Medicine & Biology 41.11 
%   (1996): 2271."
%
% INPUTS 
%       f (1_by_N double): freuency range (Hz)
%       tissue (1_by_N char): tissue name
% 
% OUTPUTS
%   gamma (1_by_N vector): complex propagation constant of medium
%   sigma (1_by_N vector): total conductivity of medium (S/m)
%   eps (1_by_N vector): real part of mediums' complex permittivity (F/m)

e0=8.854187817e-12;
mu0 = 4*pi*1e-7;
w = 2 * pi * f;

switch tissue
    case 'Blood'
        input = [4.0 56 8.38 0.10 5200 132.63 0.10 0 0 0 0 0 0 0.7000];
    case 'Bone_cancellous'
        input = [ 2.5 18.0 13.26 0.22 300 79.58 0.25 2e4 159.15 0.20 2e7 15.915 0.00 0.0700];
    case 'Bone_cortical'
        input = [2.5 10.0 13.26 0.20 180 79.58 0.20 5e3 159.15 0.20 1e5 15.915 0.00 0.0200];
    case 'Brain_grey_matter'
        input =  [4.0 45.0 7.96 0.10 400 15.92 0.15 2e5 106.10 0.22 4.5e7  5.305 0.00 0.0200];
    case 'Brain_white_matter'
        input =  [4.0 32.0 7.96 0.10 100 7.96 0.10 4e4 53.05 0.30 3.5e7 7.958 0.02 0.0200];
    case 'Fat_infiltrated'
        input =  [2.5 9.0 7.96 0.20 35 15.92 0.10 3.3e4 159.15 0.05 1e7 15.915 0.01 0.0350];
    case 'Fat_not_infiltrated'
        input = [2.5 3.0 7.96 0.20 15 15.92 0.10 3.3e4 159.15 0.05 1e7 7.958 0.01 0.0100];
    case 'Heart'
        input =  [4.0 50.0 7.96 0.10 1200 159.15 0.05 4.5e5 72.34 0.22 2.5e7 4.547 0.00 0.0500];
    case 'Kidney'
        input =  [4.0 47.0 7.96 0.10 3500 198.94 0.22 2.5e5 79.58 0.22 3e7 4.547 0.00 0.0500];
    case 'Lens_cortex'
        input =  [4.0 42.0 7.96 0.10 1500 79.58 0.10 2e5 159.15 0.10 4e7 15.915 0.00 0.3000];
    case 'Liver'
        input =  [4.0 39.0 8.84 0.10 6000 530.52 0.20 5e4 22.74 0.20 3e7 15.915 0.05 0.0200];
    case 'Lung_inflated'
        input = [ 2.5 18.0 7.96 0.10 500 63.66 0.10 2.5e5 159.15 0.20 4e7  7.958 0.00 0.0300];
    case 'Muscle'
        input = [4.0 50 7.23 0.10 7000 353.68 0.10 1.2e6 318.31 0.10 2.5e7 2.274 0 0.2];
    case 'Skin_dry'
        input =  [ 4.0 32.0 7.23 0.00 1100 32.48 0.20 0 0 0 0 0 0  0.0002 ];
    case 'Skin_wet'
        input =  [4.0 39.0 7.96 0.10 280 79.58 0.00 3e4 1.59 0.16 3e4 1.592 0.20 0.0004];
    case 'Spleen'
        input =  [4.0 48.0 7.96 0.10 2500 63.66 0.15 2e5 265.26 0.25 5e7 6.366 0.00 0.0300];
    case 'Tendon'
        input = [4.0 42.0 12.24 0.10 60 6.37 0.10 6e4 318.31 0.22 2e7 1.326 0.00 0.2500];
end

e_inf = input(1,1);
delta_e = input(1,2:3:11);
tau = [input(1,3)*1e-12 input(1,6)*1e-9 input(1,9)*1e-6 input(1,12)*1e-3];
alpha = input(1,4:3:13);
sigmai = input(1,14);

% Cole-Cole Approximation
e = e_inf + sigmai./(1i*w*e0);
for n = 1:4
    e = e + delta_e(n)./( 1+(1i*w*tau(n)).^(1-alpha(n)));
end

eps = real(e)*e0;
sigma = -imag(e) .* w * e0 ;
gamma = sqrt( 1i*2*pi.*f*mu0.*(sigma + 1i*2*pi*f.*eps) );

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%