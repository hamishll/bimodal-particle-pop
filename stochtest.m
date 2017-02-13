clear all; close all;

%% Generating particles
% Args: N, Mean size, Variance, P(center), P(edge)
[ParticleA, InputsA] = StochasticS(200, 50, 5, 1, 0.1);
[ParticleB, InputsB] = StochasticS(1000, 10, 2, 0.1, 1);

%% Plotting data
figure(1); set(figure(1), 'Position', [30, 70, 900, 900]);

subplot(2,2,1); title('Bimodal flake size distributions');
hold on
    x = [0:1:100];
    normA = normpdf(x,InputsA(2),InputsA(3));
    normB = normpdf(x,InputsB(2), InputsB(3));
    plot(x,normA, 'r');
    plot(x,normB , 'b');
hold off

subplot(2,2,2); title('Linear distribution');
hold on
    xl = [0:0.05:1];
    distA = (InputsA(5)-InputsA(4))*xl+InputsA(4);
    distB = (InputsB(5)-InputsB(4))*xl+InputsB(4);
    plot(xl, distA, 'r');
    plot(xl, distB, 'b');
    axis([0,inf,0,inf]);
hold off

subplot(2,2,3); title('Radial distribution');
hold on
    xl = [0:0.05:1];
    distA = (InputsA(5)-InputsA(4))*(xl.^2)+xl*InputsA(4);
    distB = (InputsB(5)-InputsB(4))*(xl.^2)+xl*InputsB(4);
    plot(xl, distA, 'r');
    plot(xl, distB, 'b');
    axis([0,inf,0,inf]);
    scatter([ParticleA.x], [ParticleA.y], [ParticleA.s], 'r');
    scatter([ParticleB.x], [ParticleB.y], [ParticleB.s], 'b');
hold off

subplot(2,2,4); title('2D distribution');
hold on
    scatter([ParticleA.X], [ParticleA.Y], [ParticleA.s], 'r');
    scatter([ParticleB.X], [ParticleB.Y], [ParticleB.s], 'b');
hold off