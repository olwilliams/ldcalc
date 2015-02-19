function [results] = ldcalc(parameters)
%Calculates linear dichroism for a patch on a membrane. Parameter options:
%sphereposition, probe1, probe2, excite, dipole_angle, number_particles,
%number_steps, stepscale, showplots.

disp('Running, please wait...')
%Default parameters
default.sphereposition = [1;0;90];
default.probe1 = [0;0;1]; %Parallel - z direction
default.probe2 = [1;0;0]; %Perpendicular - x direction
default.excite = [0;0;1]; %Vertical excitation
default.dipole_angle = 45;
default.number_particles = 100;
default.number_steps = 15000;
default.stepscale = 1;
default.showplots = 1;

%Check for existance of input parameters; use default if none found
if exist('parameters','var') == 0
    parameters = default;
end
if isfield(parameters,'sphereposition') == 0
    parameters.sphereposition = default.sphereposition;
end
if isfield(parameters,'probe1') == 0
    parameters.probe1 = default.probe1;
end
if isfield(parameters,'probe2') == 0
    parameters.probe2 = default.probe2;
end
if isfield(parameters,'excite') == 0
    parameters.excite = default.excite;
end
if isfield(parameters,'dipole_angle') == 0
    parameters.dipole_angle = default.dipole_angle;
end
if isfield(parameters,'number_particles') == 0
    parameters.number_particles = default.number_particles;
end
if isfield(parameters,'number_steps') == 0
    parameters.number_steps = default.number_steps;
end
if isfield(parameters,'stepscale') == 0
    parameters.stepscale = default.stepscale;
end
if isfield(parameters,'showplots') == 0
    parameters.showplots = default.showplots;
end

%Invoke patchdiffuse to calculate LD
sphereposition = parameters.sphereposition;
probe1 = parameters.probe1;
probe2 = parameters.probe2;
excite = parameters.excite;
dipole_angle = parameters.dipole_angle;
number_particles = parameters.number_particles;
number_steps = parameters.number_steps;
stepscale = parameters.stepscale;
showplots = parameters.showplots;

[parprobe, perpprobe] = patchdiffuse(sphereposition,probe1,probe2,excite,dipole_angle,number_particles,number_steps,stepscale,showplots);

%Compile results of simulation
LD = parprobe - perpprobe;
isotropic = (parprobe + 2*perpprobe)/3;
anisotropy = LD./isotropic;
[tau_1exp,amp_1exp,constant_1exp,model_1exp,exitflag1exp] = onetaufit(LD);
[tau_2exp,amp_2exp,constant_2exp,model_2exp,exitflag2exp] = twotaufit(LD);
results.LD = LD;
results.anisotropy = anisotropy;
results.tau1exp = tau_1exp;
results.amp1exp = amp_1exp;
results.constant1exp = constant_1exp;
results.model1exp = model_1exp;
results.tau2exp = tau_2exp;
results.amp2exp = amp_2exp;
results.constant2exp = constant_2exp;
results.model2exp = model_2exp;

%Display results if showplots is true
if showplots == 1;
    figure
    hold on
    plot(LD,'b')
    plot(anisotropy,'g')
    xlabel('Number of steps')
    title('Simulation')
    legend('Linear Dichroism','Anisotropy')
    
    figure
    hold on
    plot(LD,'b')
    plot(model_1exp,'g')
    plot(model_2exp,'r')
    xlabel('Number of steps')
    ylabel('LD')
    title('Simulated LD and exponential models')
    legend('LD','1 exp fit','2 exp fit')
end

fprintf('\n~~~~~~~~Model parameters~~~~~~~~\n')
fprintf('Number of dipoles: %d\nNumber of steps: %d\nStepping angle: %f\nDipole_angle: %f\nPosition on sphere: r = %1.3f, theta = %1.3f, phi = %1.3f\n',number_particles,number_steps,stepscale,dipole_angle,sphereposition(1),sphereposition(2),sphereposition(3))
fprintf('\n~~~~~~~~Fitting Results~~~~~~~~\n')
if exitflag1exp == 0
    fprintf(2,'\nWarning: Single exponential fitting did not converge!')
end
fprintf('\nSingle exponential:\nRate constant: %f\nAmplitude: %f\nConstant term: %f\n',tau_1exp,amp_1exp,constant_1exp)
if exitflag2exp == 0
    fprintf(2,'\nWarning: Double exponential fitting did not converge!')
end
fprintf('\nDouble exponential:\nRate constant 1: %f\nAmplitude 1: %f\nRate constant 2: %f\nAmplitude 2: %f\nConstant term: %f\n',tau_2exp(1),amp_2exp(1),tau_2exp(2),amp_2exp(2),constant_2exp)
end

function [probe1, probe2] = patchdiffuse(sphereposition,probe1dir,probe2dir,excitedir,dipole_angle,number_particles,number_steps,stepscale,showplots)
%Rotational diffusion function. Starts with a location on a unit sphere,
%sets up a molecular dipole at an angle to the plane normal, and picks the
%part of the cone with maximal overlab with a vertical excitation. This
%initial vector is then perturbed by a normally distributed angle about the
%plane normal via a rotation matrix.
%Convention I'm using: (r,theta,phi) where theta is the azimuthal angle in
%the xy plane [0,2pi) and phi is the polar angle [0,pi].

%Initialization
rng('shuffle')
absorbance = zeros(2,number_steps);

%Make sure directions are unit vectors
probe1unit = probe1dir/norm(probe1dir);
probe2unit = probe2dir/norm(probe2dir);
exciteunit = excitedir/norm(excitedir);

%Convert sphere location to Cartesian unit vector
planenormal = spheretocart([1;sphereposition(2:3)]);

%Initialize dipole matrix
%Define dipole location; the dipole angle is added and subtracted from the
%membrane normal vector. I pick the one that has the most projection along
%the excitation direction (probably the z direction).
dipolevectorpositive = [sphereposition(1);sphereposition(2);sphereposition(3)+dipole_angle];
dipolevectorpositive(2) = validatetheta(dipolevectorpositive(2));
dipolevectorpositive(3) = validatephi(dipolevectorpositive(3));
dipolevectorpositive_cart = spheretocart(dipolevectorpositive);
dipolevectornegative = [sphereposition(1);sphereposition(2);sphereposition(3)-dipole_angle];
dipolevectornegative(2) = validatetheta(dipolevectornegative(2));
dipolevectornegative(3) = validatephi(dipolevectornegative(3));
dipolevectornegative_cart = spheretocart(dipolevectornegative);
overlap_pos = (dot(dipolevectorpositive_cart,exciteunit))^2;
overlap_neg = (dot(dipolevectornegative_cart,exciteunit))^2;
if overlap_pos >= overlap_neg
    dipolevector = dipolevectorpositive_cart;
elseif overlap_neg > overlap_pos
    dipolevector = dipolevectornegative_cart;
end

%Create dipole matrix
dipolematrix = repmat(dipolevector,1,number_particles);

%Diffuse the dipoles
for timestep = 1:1:number_steps
    for currentdipole = 1:1:number_particles
        perturbation = randn*stepscale;
        rotatematrix = rotateaxis(planenormal,perturbation);
        dipolematrix(:,currentdipole) = rotatematrix*dipolematrix(:,currentdipole);
    end
    %Capture absorbance changes per timepoint
    absorbance(1,timestep) = sum(sum((repmat(probe1unit,1,number_particles).*dipolematrix)).^2);
    absorbance(2,timestep) = sum(sum((repmat(probe2unit,1,number_particles).*dipolematrix)).^2);
end

%Compute results
probe1 = absorbance(1,:);
probe2 = absorbance(2,:);

%Display distribution if showplots is true
if showplots == 1
    figure
    hold on
    for n = 1:number_particles
        xline = dipolematrix(1,n);
        yline = dipolematrix(2,n);
        zline = dipolematrix(3,n);
        plot3([0 xline],[0 yline],[0 zline])
        plot3([0 -xline],[0 -yline],[0 -zline])
    end
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view([0 0])
    title('Final Dipole Distribution')
end
end

function [newtheta] = validatetheta(oldtheta)
%Make sure theta (in spherical coordinates) is in the interval [0,360)
%degrees.
theta = oldtheta;
while theta >= 360
    theta = theta - 360;
end
while theta < 0
    theta = theta + 360;
end
newtheta = theta;
end

function [newphi] = validatephi(oldphi)
%Check that phi is in the interval [0,180] degrees.
phi = oldphi;
while phi > 180
    phi = phi - 360;
end
while phi < -180
    phi = phi + 360;
end
if (phi < 0)&&(phi > -180)
    phi = abs(phi);
end
newphi = phi;
end
        
function [cartesianvector] = spheretocart(sphericalvector)
%Subfunction to convert a spherical coordinate vector to Cartesian. 
%Convention I'm using: (r,theta,phi) where theta is the azimuthal angle in
%the xy plane [0,2pi) and phi is the polar angle [0,pi). Theta is
%originally in degrees but is converted to radians.
r = sphericalvector(1);
theta = sphericalvector(2)*(pi/180);
phi = sphericalvector(3)*(pi/180);
x = r*cos(theta)*sin(phi);
y = r*sin(theta)*sin(phi);
z = r*cos(phi);
cartesianvector = [x;y;z];
end

function [spherevector] = carttosphere(cartesianvector)
%Subfunction to convert Cartesian vectors to spherical coordinates. Angular
%values are converted into degrees.
x = cartesianvector(1);
y = cartesianvector(2);
z = cartesianvector(3);
r = sqrt(x^2+y^2+z^2);
theta = atan2(y/x);
if theta < 0
    theta = theta + 2*pi;
end
theta = theta*(180/pi);
phi = acos(z/r);
phi = phi*(180/pi);
spherevector = [r;theta;phi];
end

function [rotationmatrix] = rotateaxis(directionvector,degangle)
%Subfunction to calculate the rotation matrix about an arbitrary axis. The
%direction vector should be in Cartesian coordinates with length 1. The
%angle is input in degrees but is then converted to radians to take
%advantage of cos() being significantly faster than cosd().

ux = directionvector(1);
uy = directionvector(2);
uz = directionvector(3);
radangle = degangle*(pi/180);

R11 = cos(radangle) + ux^2*(1-cos(radangle));
R12 = ux*uy*(1-cos(radangle)) - uz*sin(radangle);
R13 = ux*uz*(1-cos(radangle)) + uy*sin(radangle);
R21 = uy*ux*(1-cos(radangle)) + uz*sin(radangle);
R22 = cos(radangle) + uy^2*(1-cos(radangle));
R23 = uy*uz*(1-cos(radangle)) - ux*sin(radangle);
R31 = uz*ux*(1-cos(radangle)) - uy*sin(radangle);
R32 = uz*uy*(1-cos(radangle)) + ux*sin(radangle);
R33 = cos(radangle) + uz^2*(1-cos(radangle));
rotationmatrix = [R11 R12 R13;R21 R22 R23;R31 R32 R33];
end

function [tau,amp,constant,model,exitflag] = onetaufit(lddata)
%Subfunction to fit a single exponential curve
number_timesteps = length(lddata);
times = 1:1:number_timesteps;
tauguess = number_timesteps/2;
ampguess = lddata(round(number_timesteps/2));
constantguess = lddata(end);
guesses = [tauguess;ampguess;constantguess];
[results,~,exitflag] = fminsearch(@(x) onetaumin(x,lddata),guesses);
tau = results(1);
amp = results(2);
constant = results(3);
model = amp*exp(-times./tau)+constant;
end

function [residual] = onetaumin(guessparameters,data)
%Minimization functional for onetaufit
number_timesteps = length(data);
times = 1:1:number_timesteps;
tau = guessparameters(1);
amp = guessparameters(2);
constant = guessparameters(3);
modeldata = amp*exp(-times./tau) + constant;
difference = data - modeldata;
squaredifference = difference.^2;
residual = sum(squaredifference);
end

function [tau,amp,constant,model,exitflag] = twotaufit(lddata)
%Subfunction to fit a double exponential curve
number_timesteps = length(lddata);
times = 1:1:number_timesteps;
tauguess1 = number_timesteps*(1/3);
tauguess2 = number_timesteps*(2/3);
ampguess1 = lddata(round(number_timesteps/3));
ampguess2 = lddata(round(number_timesteps*2/3));
constantguess = lddata(end);
guesses = [tauguess1;tauguess2;ampguess1;ampguess2;constantguess];
[results,~,exitflag] = fminsearch(@(x) twotaumin(x,lddata),guesses);
tau = results(1:2);
amp = results(3:4);
constant = results(5);
model = amp(1)*exp(-times./tau(1)) + amp(2)*exp(-times./tau(2)) + constant;
end

function [residual] = twotaumin(guessparameters,data)
%Minimization functional for onetaufit
number_timesteps = length(data);
times = 1:1:number_timesteps;
tau1 = guessparameters(1);
tau2 = guessparameters(2);
amp1 = guessparameters(3);
amp2 = guessparameters(4);
constant = guessparameters(5);
modeldata = amp1*exp(-times./tau1) + amp2*exp(-times./tau2) + constant;
difference = data - modeldata;
squaredifference = difference.^2;
residual = sum(squaredifference);
end
