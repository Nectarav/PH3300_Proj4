mu0 = 4*pi*1e-7;
% Constant mu0 

I = 5; 
R = 0.1;
Z = linspace(0, 0.5, 100);

N = [20, 50, 100, 500];

Bval = zeros(numel(N), numel(Z), 3);
% Stores numerical values of the magnetic field for N, z, and 
% its vector coordinates respectively
Bactualval = zeros(numel(N), numel(Z), 3);
% Stores actual values of the magnetic field for N, z, and 
% its vector coordinates respectively
Berror = zeros(numel(N), numel(Z), 3);
% Same, but for the difference between B numerical and actual
Errorrel = zeros(numel(N), numel(Z));
% Values for relative error for N, z. Note that this is a scalar quantity  


for k = 1:length(N)
    n = N(k);
    dtheta = (2*pi)/n;

    for j = 1:length(Z)
        z = Z(j);
        B = [0,0,0];

        for i = 0:n-1
            theta1 = i * dtheta;
            theta2 = (i+1) * dtheta;
            r = [-R * cos(theta1), -R * sin(theta1), z];
            l1 = [R * cos(theta1), R * sin(theta1), 0];
            l2 = [R * cos(theta2), R * sin(theta2), 0];
            dl = (l2 - l1);
            dBi = ((mu0 * I/(4 * pi)) * cross(dl, r) / norm(r)^3); 
            B = B + dBi;
            % Calculating dBi and adding it to the total magnetic field
        end 

        Bzactual = (mu0 * I * R^2)/(2*(R^2 + z^2)^(3/2)) * [0,0,1]; 
        % Computing the actual magnetic field based on our derivation


        Bval(k, j, :) = B;
        Bactualval(k, j, :) = Bzactual;
        Berror(k, j, :) = B - Bzactual;
        Errorrel(k, j) = norm(B - Bzactual) / norm(Bzactual) * 100;
    end 
    figure
    plot(Z, squeeze(Bval(k,:,3)), 'blue', 'DisplayName', 'Numerical');
    hold on
    plot(Z, squeeze(Bactualval(k,:,3)), 'red--', 'DisplayName', 'Actual');
    xlabel('z')
    ylabel('B_z')
    title(sprintf('B_z vs z for N = %d', n))
    legend('Location', 'best')
    grid on

    figure
    plot(Z, squeeze(Errorrel(k,:)), 'green', 'DisplayName', ...
        'Relative Error');
    xlabel('z')
    ylabel('Relative error (%)')
    ytickformat('%.2g')
    title(sprintf('Relative Error vs z for N = %d', n))
    legend('Location', 'best')
    grid on

end 

Maxerror = zeros(size(N));
for k = 1:length(N)
    for z = 1:length(Z)
        Maxerror(k) = max(Maxerror(k), norm(squeeze(Berror(k, z, :))));
    end 
end

figure
plot(N, Maxerror(:), 'blue', 'DisplayName', 'Maximum Error')
xlabel('Number of segments, N')
ylabel('Maximum absolute error')
title('Maximum Magnetic Field Error vs N')
legend('Location','best')
grid on
