function [k, info] = sader_k_from_fQ(fR, Q, varargin)
%SADER_K_FROM_FQ  Spring constant from fR & Q using Sader (rectangular beam, fluid)
%
%   k = SADER_K_FROM_FQ(fR, Q, 'L', L, 'w', w, 'rho_fl', rho, 'eta', eta, ...)
%
% Required:
%   fR        : resonance frequency [Hz]   (scalar)
%   Q         : quality factor     [-]     (scalar)
%
% Name-value (no defaults for geometry → you must supply L and w):
%   'L'       : cantilever length [m]      (e.g. 240e-6)
%   'w'       : cantilever width  [m]      (e.g.  40e-6)
%   'rho_fl'  : fluid density     [kg/m^3] (default 997 for water @ ~25°C)
%   'eta'     : dynamic viscosity  [Pa·s]  (default 0.89e-3 for water @ ~25°C)
%
% Returns:
%   k         : spring constant [N/m]
%   info      : struct with fields:
%               .omega, .Re, .Gamma, .Gamma_r, .Gamma_i, .rho_fl, .eta, .L, .w
%
% Formula (Sader liquid; rectangular beam):
%   k = 0.1906 * rho_fl * w^2 * L * Q * Gamma_i(Re) * (2*pi*fR)^2
%
% Notes:
% - This avoids thickness/material parameters by using Q and the imag part
%   of the hydrodynamic function Gamma_i(Re).

    % ------------------ Parse inputs ------------------
    p = inputParser;
    p.addRequired('fR', @(x)isscalar(x)&&isnumeric(x)&&x>0);
    p.addRequired('Q',  @(x)isscalar(x)&&isnumeric(x)&&x>0);
    p.addParameter('L',      [], @(x)isscalar(x)&&isnumeric(x)&&x>0);
    p.addParameter('w',      [], @(x)isscalar(x)&&isnumeric(x)&&x>0);
    p.addParameter('rho_fl', 997, @(x)isscalar(x)&&isnumeric(x)&&x>0);       % water
    p.addParameter('eta',  0.89e-3, @(x)isscalar(x)&&isnumeric(x)&&x>0);     % water
    p.parse(fR, Q, varargin{:});
    L      = p.Results.L;
    w      = p.Results.w;
    rho_fl = p.Results.rho_fl;
    eta    = p.Results.eta;

    if isempty(L) || isempty(w)
        error('Please provide both ''L'' (length) and ''w'' (width) in meters.');
    end

    % ------------------ Core computation ------------------
    omega = 2*pi*fR;

    % Reynolds number based on width (standard Sader definition)
    Re = (rho_fl * omega * w^2) / (4 * eta);

    % Hydrodynamic function Gamma(Re) = Omega(Re) * Gamma_circ(Re)
    [Gamma_r, Gamma_i, Gamma] = local_hydro_Gamma_rect(Re);

    % Sader liquid stiffness (rectangular beam; uses Q and Gamma_i)
    k = 0.1906 * rho_fl * (w^2) * L * Q * Gamma_i * (omega^2);

    % ------------------ Book-keeping ------------------
    if nargout > 1
        info = struct('omega', omega, 'Re', Re, ...
                      'Gamma', Gamma, 'Gamma_r', Gamma_r, 'Gamma_i', Gamma_i, ...
                      'rho_fl', rho_fl, 'eta', eta, 'L', L, 'w', w);
    end
end

% ===== Local: hydrodynamic function for rectangular beam  =====
function [Gamma_r, Gamma_i, Gamma] = local_hydro_Gamma_rect(Re)
    % Sader-style fit via Omega(Re) times circular-cylinder correction Gamma_circ(Re)

    tau = log10(Re);

    % Omega(Re): rational polynomial fits in log10(Re)
    ReOmega = polyval([ 0.00069085  -0.0035117   0.044055  -0.12886   0.46842  -0.48274   0.91324], tau) ./ ...
              polyval([ 0.00069085  -0.0035862   0.045155  -0.13444   0.48690  -0.56964   1.00000], tau);

    % Note: the first two terms are summed intentionally
    % producing a single first coefficient ~2.0067e-05 when evaluated.
    ImOmega = polyval([(-0.000044510 + 0.000064577)  -0.00010961   0.016294  -0.029256  -0.024134], tau) ./ ...
              polyval([ 0.00286361  -0.014369   0.079156  -0.18357   0.55182  -0.59702   1.00000], tau);

    Omega = ReOmega + 1i*ImOmega;

    % Circular-cylinder correction (Bessel K functions), evaluated at complex argument
    z = sqrt(1i*Re);
    Gamma_circ = 1 + (4i)*besselk(1, -1i*z) ./ ( z .* besselk(0, -1i*z) );

    Gamma   = Omega .* Gamma_circ;
    Gamma_r = real(Gamma);
    Gamma_i = imag(Gamma);
end
