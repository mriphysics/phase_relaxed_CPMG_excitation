% Function that produces "shell" trajectory proposed in Wong & Roos 
% (10.1002/mrm.1910320614).
% Produces a single shell with unit radius, M points and
% (optionally) rotated by angle rot (radians) from [1 1 1]
% function k = gen_single_shell(M,rot)
%
% Shaihan Malik 2012

function k = gen_single_shell(M,varargin)

rotate_kspace=false;
slow_poles=false;

for ii=1:length(varargin)
    % rotation
    if strcmpi(varargin{ii},'rot')
        rotate_kspace=true;
        rot = varargin{ii+1};
    end
    
    % slow down at poles
    if strcmpi(varargin{ii},'slow')
        slow_poles=true;
        alpha = varargin{ii+1};
    end
    
end

n=1:M;
c=sqrt(M*pi);

if ~slow_poles
    % Standard version
    z = (2*n-M-1)/M;
else
    % This is mentioned in the discussion of the paper (Eq.14)
    z = atan(alpha*(2*n-M-1)/M) ./ atan(alpha*(M-1)/M);
end
x = cos(c*asin(z)).*sqrt(1-z.^2);
y = sin(c*asin(z)).*sqrt(1-z.^2);

k = [x.' y.' z.'];

if rotate_kspace
    Ro = rot3d(rot*[1 1 0]/sqrt(2));
    k = k*Ro;
end

end

function R = rot3d(u)
%%% 3D rotation around vector u. u specifies rotation axis, modulus of u
%%% determines magnitude of rotation

% check zero input
if any(u)
    theta=norm(u);
    u=u/theta;
    ct=cos(theta);
    st=sin(theta);
    R = [[ct + u(1)^2*(1-ct) u(1)*u(2)*(1-ct)-u(3)*st u(1)*u(3)*(1-ct)+u(2)*st];...
        [u(2)*u(1)*(1-ct)+u(3)*st ct+u(2)^2*(1-ct) u(2)*u(3)*(1-ct)-u(1)*st];
        [u(3)*u(1)*(1-ct)-u(2)*st u(2)*u(3)*(1-ct)+u(1)*st] ct+u(3)^2*(1-ct)];

else
    R=eye(3);
end

end