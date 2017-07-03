%  Input
% initialization - trajectory_generator([], [], waypoints)
% runtime - trajectory_generator(t, state, []) 

% ---- State
%       state: The current state of the robot with the following fields:
%       state.pos = [x; y; z], 
%       state.vel = [x_dot; y_dot; z_dot],
%       state.rot = [phi; theta; psi], 
%       state.omega = [p; q; r]

% ---- t
%       time


%  Output
%       des_state.pos = [x; y; z], 
%       des_state.vel = [x_dot; y_dot; z_dot],
%       des_state.acc = [x_ddot; y_ddot; z_ddot], 
%       des_state.yaw,
%       des_state.yawdot


% Matrix schemas
%
%  ----- s_derivative_matrix
%           R t^7, t^6, ... t^0 
%           V t^7, t^6, ... t^0
%           A t^7, t^6, ... t^0
%           J t^7, t^6, ... t^0
%
%  ----- alphaCoefficients
%           rInit   t^7, t^6, ... t^0 
%           rFinal  t^7, t^6, ... t^0 
%           vInit   t^7, t^6, ... t^0 
%           vFinal  t^7, t^6, ... t^0 
%           aInit   t^7, t^6, ... t^0 
%           aFinal  t^7, t^6, ... t^0 
%           jInit   t^7, t^6, ... t^0 
%           jFinal  t^7, t^6, ... t^0 


%waypoints = [0,   0,   0;
%             1,   1,   1;
%             2,   0,  2;
%             3,   -1,  1;
%             4,   0,   0]
function [ desired_state ] = traj_generator(t, state, waypoints)


persistent waypoints0 traj_time d0 pos_coeffs vel_coeffs acc_coeffs

if nargin > 2
    number_waypoints = size(waypoints);
number_waypoints = number_waypoints(1)
number_functions = number_waypoints - 1

% Create symbolic velocities, acc, jerks
midpoint_velocities = []

tP = 5;
tMax = tP * number_functions;

for AXIS = 1:1
    
    % Matrices
    s_derivative_matrix = getSDerivativeMatrix()
    alphaCoefficients = getAlphaCoefficients(number_functions,s_derivative_matrix)
    values = getValues( waypoints, AXIS,number_functions)
    
    alphas = linsolve(alphaCoefficients,values)

    derivative = 0;
    print_general_equation = false;

end
          

else
       des_state.pos = [x(t); y(t); z(t)]; 
       des_state.vel = [x_dot(t); y_dot(t); z_dot(t)];
       des_state.acc = [x_ddot(t); y_ddot(t); z_ddot(t)];
       des_state = 0;
       des_state = 0;
end

end


function [values] = getValues(waypoints,axis,number_functions)
    
    values = sym(zeros(number_functions*8,1));
    
    for function_inc = 0:number_functions-1
        
        % r
        values( (function_inc * 8) + 1, 1) =  waypoints(function_inc+1,axis);        
        values((function_inc * 8) + 2 , 1) = waypoints(function_inc+2,axis);
        
        % v
        values( function_inc * 8 + 3, 1) = midpoint_velocities(function_inc+1,axis);
        values( function_inc * 8 + 4, 1) = midpoint_velocities(function_inc+2,axis);

        % a
        values( function_inc * 8 + 5, 1) = midpoint_accellerations(function_inc+1,axis);
        values( function_inc * 8 + 6, 1) = midpoint_accellerations(function_inc+2,axis);

        % j
        values( function_inc * 8 + 7, 1) = midpoint_jerks(function_inc+1,axis);
        values( function_inc * 8 + 8, 1) = midpoint_jerks(function_inc+2,axis);
    end
end

function [alphaCoefficients] = getAlphaCoefficients(number_functions,s_derivative_matrix)
    % Set alphaCoefficients
    alphaCoefficients = sym( zeros( number_functions * 8, 8 ));
    for function_num = 0:number_functions-1
    
        for row = 1:4

            for col = 1:8
                
                % alpha coefficient row
                alpha_row = (function_num * 8 ) + (row*2) - 1;
                
                % set init
                alphaCoefficients(alpha_row,col) = s_derivative_matrix(row,col);
                
                % set final
                alphaCoefficients(alpha_row+1,col) = s_derivative_matrix(row,col);

            end
        end
    end

end





% derivatives
%   0 - position
%   1 - vel
%      ...
function [s] = getS(alphas,t,Tp,Ts,derivative,print_general_equation)
    s = 0;
    
    Tpiece = Tp;
    Tstart_i = Ts;
    
    if print_general_equation
        fprintf('p%1.0f(t):',derivative)
    end
    
    for row = 1:8
        
        t_coeff = subs(alphas(row,1));
        
        syms t_temp
        s_term = ( (t_temp - Ts)/(Tp) ) ^ (8-row);
        
        if derivative>0
            for i = 0:derivative-1
                s_term = sym(simplify(diff( s_term,t_temp)));
            end
        end
        
        s_coeff = double(subs(s_term,t_temp,1));
        s_term = subs(s_term,t_temp,t);
        
        if abs(t_coeff) > 0
            
            if print_general_equation
                fprintf( ' (%6.2f * (%6.10f * t^%1.0f))+', t_coeff,s_coeff,(8-row-derivative));
            end
            
            s = s + t_coeff * s_term;
        end
    end
    
    if print_general_equation
        fprintf('\n\n')
    end
    
    s = simplify(s);
end




function [s_derivative_matrix] = getSDerivativeMatrix()
    
    syms s t Tstart_i Tpiece
    s = (t - Tstart_i)/Tpiece;
    
    s_derivative_matrix = sym(zeros(4,8));
    
    s_derivative_matrix(1,8) = 1;
    
    % Set first row (position)
    for col = 1:7
        s_derivative_matrix(1,col) = sym(   s^(8-col)  );
    end
    
    % Set velocity, acceleration, jerk rows
    for row = 2:4
        for col = 1:8
            s_derivative_matrix(row,col) = sym(simplify(diff(s_derivative_matrix(row-1,col),t)));
        end    
    end
end



        
        