% Title: 3D visualization of Grover's Search
% Description: 
% Visualizing operations on sphere 
% Author: Amritesh Sharma

clear all;
close all;

% draw the circle
hold on
draw_sphere();
draw_circle(0,0,0);

%% desired state
theta_desired = pi/2;
phi_desired = pi/2;
r1 = sin(phi_desired);
desired_state_pt = [sin(theta_desired)*cos(phi_desired) sin(theta_desired)*sin(phi_desired) cos(theta_desired)];

%% initial state
theta_init = pi/2;
phi_init = 0.3*pi/2;
r_init = cos(phi_init);

origin_pt = [0 0 0];                         % First Point
init_state_pt = [sin(theta_init)*cos(phi_init) sin(theta_init)*sin(phi_init) cos(theta_init)];   % Second Point
dp = init_state_pt-origin_pt;                         % Difference
quiver3(origin_pt(1),origin_pt(2),origin_pt(3),dp(1),dp(2),dp(3),0,'g','LineWidth',2)

%% final stage

theta_final = pi/2;
phi_final = phi_init;
r1 = sin(phi_final);
% draw the circle
draw_circle(r1,theta_final,phi_final);

%% Iterations
N = floor((pi/2-phi_init)/(2*phi_init));
phi_iter = phi_init;
Path = [];
path_start_idx = 1;

for j = 1:N

% step 1: oracle
Path = draw_circle(cos(phi_iter), theta_init, 0, init_state_pt', desired_state_pt', Path,'oracle');
path_start_idx = draw_path(Path, path_start_idx);

% step 2: inversion about mean
Path = draw_circle(cos((phi_iter+2*phi_init)-phi_init), theta_init, phi_init, init_state_pt', desired_state_pt', Path, 'inversion');
path_start_idx = draw_path(Path, path_start_idx);

phi_iter = (phi_iter+2*phi_init);
p2 = [sin(theta_init)*cos(phi_iter) sin(theta_init)*sin(phi_iter) cos(theta_init)];   % Second Point
dp = p2-origin_pt;                         % Difference
quiver3(origin_pt(1),origin_pt(2),origin_pt(3),dp(1),dp(2),dp(3),0,'b--','LineWidth',3)

end

% step 1: oracle
Path = draw_circle(cos(phi_iter), theta_init, 0, init_state_pt', desired_state_pt', Path,'oracle');
path_start_idx = draw_path(Path, path_start_idx);

final_intersection_coordinate = inv(YawPitchMatrix(phi_init,theta_init))*Path(:,end);
err_angle = -abs(acos(final_intersection_coordinate(2)));

% step 2: inversion about mean
Path = draw_circle(r1, theta_init, phi_init, init_state_pt', desired_state_pt', Path, 'inversion', err_angle);
path_start_idx = draw_path(Path, path_start_idx);

%% Custom defined functions
function y = draw_sphere()
    hold on;
    [X,Y,Z] = sphere(100);
    surf(X,Y,Z,'FaceAlpha',0.5,'EdgeAlpha',0.1, 'EdgeColor','k')
    view(120,30)
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis equal

    % draw the x-axis
    p1 = [-1.5 0 0];                         % First Point
    p2 = [1.5 0 0];                         % Second Point
    dp = p2-p1;                         % Difference
    quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0,'k','LineWidth',2)

    % draw the y-axis
    p1 = [0 -1.5 0];                         % First Point
    p2 = [0 1.5 0];                         % Second Point
    dp = p2-p1;                         % Difference
    quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0,'k','LineWidth',2)

    % draw the z-axis
    p1 = [0 0 -1.5];                         % First Point
    p2 = [0 0 1.5];                         % Second Point
    dp = p2-p1;                         % Difference
    quiver3(p1(1),p1(2),p1(3),dp(1),dp(2),dp(3),0,'k','LineWidth',2)

    grid minor
end

function y = draw_circle(r,theta, phi, init_state, desired_state, Path, flag, err_angle)

    % r = distance of the plane intersecting sphere
    % theta = polar angle of r
    % phi = azimuthal angle
    % flag = 'inversion' for inversion about a state operation
    % flag = 'oracle' for oracle operation
    % theta_init = polar angle of initial state
    % phi_init = azimuthal angle of initial state
    switch nargin
        case 3
            flag = 'inversion';
            Path = [];
            init_state = [0 0 0]';
            desired_state = [0 0 0]';
            error_amount = -2*pi;
        case 8
            error_amount = err_angle;
        otherwise
            error_amount = -2*pi;
    end
%     if nargin <= 3
%         flag = 'inversion';
%         Path = [];
%     end

%     theta = pi/2;
%     phi = 0.2*pi/2;
%     r = sin(phi);

    pitch_angle = -(pi/2 - theta);

%%%%  Reference for drawing drawing arbitrary oriented circles - https://in.mathworks.com/matlabcentral/answers/504095-draw-a-circle-for-arbitrary-orientation-on-spherical-surface
%     YAW = [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1]; % Planar YAW rotation
%     PITCH = [cos(pitch_angle), 0, sin(pitch_angle); 0,1,0; -sin(pitch_angle), 0, cos(pitch_angle)]; % Planar PITCH rotation
%     YP = YAW*PITCH; % YAW-PITCH rotation matrix
    YP = YawPitchMatrix(phi,pitch_angle);

    r_circle = sqrt(1-r^2);
    theta_circle = [0:-0.04:-pi -pi (-pi-0.1):-0.04:-2*pi -2*pi];
    for j = 1:numel(theta_circle)
        Circle(:,j) = YP*[r r_circle*cos(theta_circle(j)) r_circle*sin(theta_circle(j))]';
        if theta_circle(j) >= -pi && flag == "oracle" 
            abs(dot(Circle(:,j),init_state))
            if abs(dot(Circle(:,j),init_state)) - init_state(2) > 1e-2
                Path = [Path Circle(:,j)];
            end
        elseif theta_circle(j) <= max(-pi,error_amount) && flag == "inversion"
            if abs(abs(dot(Circle(:,j),desired_state)) - 1) > 1e-2
                 Path = [Path Circle(:,j)];
            end
        end
    end
    
    plot3(Circle(1,:),Circle(2,:),Circle(3,:),'LineWidth',2);
    y = Path;

end

function YP = YawPitchMatrix(phi,pitch_angle)
    YAW = [cos(phi), -sin(phi), 0; sin(phi), cos(phi), 0; 0, 0, 1]; % Planar YAW rotation
    PITCH = [cos(pitch_angle), 0, sin(pitch_angle); 0,1,0; -sin(pitch_angle), 0, cos(pitch_angle)]; % Planar PITCH rotation
    YP = YAW*PITCH; % YAW-PITCH rotation matrix
end

function y = draw_path(Path, path_start_idx)

    for jj = path_start_idx : numel(Path(1,:))
        children = get(gca, 'children');
        % If there is only one, this is the first iteration, so do nothing.  Otherwise, remove the line
        if length(children)>1 && jj>path_start_idx
            delete(children(2));
        end

        line([0 Path(1,jj)],[0 Path(2,jj)],[0 Path(3,jj)], 'LineWidth', 2, 'Color', [0 1 0])
        plot3(Path(1,path_start_idx:jj),Path(2,path_start_idx:jj),Path(3,path_start_idx:jj),'r','LineWidth',3)
        pause(0.01)
        drawnow;
    end
    
    y = length(Path(1,:))-1;
end