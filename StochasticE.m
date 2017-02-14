function [ s, inputs ] = StochasticE( N , N_mean, N_var, a )
%StochasticS Summary of this function goes here
%   Detailed explanation goes here

% writing the inputs to an array
inputs = [N , N_mean, N_var, a];

% we need to multiply our input/observed probability densities
% by x in order to get to the radial dist we use for plotting

% declaring empty arrays
s = struct('x', [], 'y', [], 'a', [], 's', [], 'X', [], 'Y', []);

quadrant = 1;
count = 1;

while count < (N+1)
    
    % LINEAR
    rng_x = rand(1);
    
    % EXPONENTIAL
    y_eqn = exp(-a*(1-rng_x)) * rng_x;
    y_max = 1 * exp(0);
    
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

