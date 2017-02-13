function [ s, inputs ] = StochasticS( N , N_mean, N_var, Y1, Y2 )
%StochasticS Summary of this function goes here
%   Detailed explanation goes here

% writing the inputs to an array
inputs = [N , N_mean, N_var, Y1, Y2];

% we need to multiply our input/observed probability densities
% by x in order to get to the radial dist we use for plotting

%linear eqn values
m = Y2-Y1;
c = Y1;

% declaring empty arrays
s = struct('x', [], 'y', [], 'a', [], 's', [], 'X', [], 'Y', []);

quadrant = 1;
count = 1;

while count < (N+1)
    
    % LINEAR
    rng_x = rand(1);
    y_eqn = (m*rng_x + c)*rng_x;
    
    % EXPONENTIAL
    %a = 5;
    %y = 1 * exp(-a*(xlim-rng_x));
    %y_max = 1 * exp(0);
    
    rng_y = rand(1);
    rng_a = rand(1)*2*pi/4;
    rng_s = N_mean + randn(1)*N_var;
    
    % evaluating whether values fall within probability density
    if (rng_y < y_eqn)
        s(count).x = rng_x;
        s(count).y = rng_y; % this line is irrelevant to 2D distribution
        s(count).a = rng_a;
        s(count).s = rng_s;
        
        h = rng_x;  % h is the hypoteneuse
        dx = h - h*cos(rng_a);
        dy = h*sin(rng_a);
        
        % assigning points across each quadrant
        if quadrant == 1
            s(count).X = (h - dx);
            s(count).Y = (+dy);
        
        elseif quadrant == 2
            s(count).X = (-h + dx);
            s(count).Y = (+dy);
        
        elseif quadrant == 3
            s(count).X = (-h + dx);
            s(count).Y = (-dy);
        
        elseif quadrant == 4
            s(count).X = (+h - dx);
            s(count).Y = (-dy);
        end
        
        % increment counter
        if quadrant == 4
            quadrant = 1;
        else
            quadrant = quadrant + 1;
        end
        count = count + 1;
    end
end

end

