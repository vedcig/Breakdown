classdef Calculate
    %CALCULATE Calculate breakdown voltage
    
    properties
        d7
        p
        voltage
        radius1
        radius2
    end
    
    methods
        function [outputArg1, outputArg2] = method1(obj)
            
            obj.radius1 = obj.radius1 / 1000;
            obj.radius2 = obj.radius2 / 1000;
            V2 = 0; % V
            distance = 0.001 : 0.001 : 0.1; % m
            d = obj.radius1 + obj.radius2 + distance; % m
            eps0 = 1;
            
            N = 15;
            
            q1 = zeros(N, size(distance, 2));
            q2 = zeros(N, size(distance, 2));
            
            q1(1, :) = 4*pi*eps0*obj.radius1*1;
            q2(1, :) = 4*pi*eps0*obj.radius2*V2;
            
            b1 = zeros(N, size(distance, 2));
            b2 = zeros(N, size(distance, 2));
            
            d1 = zeros(N, size(distance, 2));
            d2 = zeros(N, size(distance, 2));
            d1(1, :) = d;
            d2(1, :) = d;

            ymin = - distance(end) * 2;
            ymax = distance(end) * 2;
            xmin = - (obj.radius1 + distance(end));
            xmax = obj.radius1 + distance(end);

            step = 0.00025; % m

            x = single(xmin:step:xmax);
            y = single(ymin:step:ymax);
            [X,Y] = meshgrid(x, y);
            
            aps1 = single(zeros(size(x, 2), size(y, 2), N, size(distance, 2)));
            aps2 = single(zeros(size(x, 2), size(y, 2), N, size(distance, 2)));
            
            for i = 2:N
                q1(i, :) = - obj.radius1./(d - b2(i - 1, :)).*q2(i - 1, :);
                q2(i, :) = - obj.radius2./(d - b1(i - 1, :)).*q1(i - 1, :);
                b1(i, :) = obj.radius1*obj.radius1./(d - b2(i - 1, :));
                b2(i, :) = obj.radius2*obj.radius2./(d - b1(i - 1, :));
                d1(i, :) = d - b1(i, :);
                d2(i, :) = d - b2(i, :);
            end
            
            Z = single(zeros(size(x, 2), size(y, 2), size(distance, 2)));
            
            for j = 1 : size(distance, 2)
            for i = 1:N
                aps1(:, :, i, j) = sqrt((X - b1(i, j)).^2 + Y.^2)';
                aps2(:, :, i, j) = sqrt((X - d2(i, j)).^2 + Y.^2)';
            
                Z(:, :, j) = Z(:, :, j) + 1/(4*pi*eps0)*(q1(i, j)./(aps1(:, :, i, j)) + (q2(i, j)./(aps2(:, :, i, j))));
            end
            end

            zmin = obj.radius1;
            zmax = obj.radius1 + distance;
         
            eta = zeros(1, size(distance, 2));
            
            d6 = distance; % m
            
            p1 = 0.05 : 0.05 : 0.3; % MPa
            P = size(p1, 2);

            Udpdvec = zeros(size(d6, 2), P);

            C = 1.6053e-10;
            EM = 2.165e7;
            A = 2873;
            
            K = 9.15;
            
            Ecr = 2.588e7;

            for j = 1 : size(distance, 2)
                z = zmin:step:zmax(j);
                [Ex, Ey] = gradient(-Z(:, :, j), step, step);
               
                E2 = zeros(size(z));
                
                E2(1, :) = sqrt(Ex(size(xmin:step:zmin, 2):size(xmin:step:zmax(j), 2), floor(size(y, 2) / 2)) .^ 2 + Ey(size(xmin:step:zmin, 2):size(xmin:step:zmax(j), 2), floor(size(y, 2) / 2)) .^ 2);
                  
                E3 = E2 * obj.voltage*1000;
                eta(j) = mean(E3)/ max(E3);
                
                xc = ones(P, 1);
                [~, xEmin] = min(E2);
            
                for i = 1 : P
                    lhs = K / (p1(i) * C);
                    rhs = zeros(xEmin - 1, 1);
                    for m = 2 : xEmin
                        rhs(m - 1) = Ecr^2 * trapz(z(1:m), (E2(1:m)).^2) / (E2(m))^2 - 2 * Ecr * EM * trapz(z(1:m), E2(1:m)) / E2(m) + EM^2 * (z(2) - zmin) * (m - 1) - A / C * (z(2) - zmin) * (m - 1);
                    end
                    temp = zeros(xEmin - 1, 1);
                    for m = 2 : xEmin
                        temp(m - 1) = lhs - rhs(m - 1);
                    end
                    [~, xc(i)] = min(abs(temp));
                    if (xc(i) > 1)    
                        Udpdvec(j, i) = Ecr * p1(i) / E2(xc(i) + 1);
            %         else
            %             U = Ecr * p1(i) / min(E2);
            %             m = size(z, 2);
            %             rhs = (U / p1(i))^2 * trapz(z(1:m), (E2(1:m)).^2) - 2 * EM * (U / p1(i)) * trapz(z(1:m), E2(1:m)) + EM^2 * udaljenost(j) - A / C * udaljenost(j);
            %             while (lhs > rhs)
            %                 U = U * 1.01;
            %                 rhs = (U / p1(i))^2 * trapz(z(1:m), (E2(1:m)).^2) - 2 * EM * (U / p1(i)) * trapz(z(1:m), E2(1:m)) + EM^2 * udaljenost(j) - A / C * udaljenost(j);
            %             end
            %             Udpdvec(j, i) = U;
                    end
                end
            end
            outputArg1 = Udpdvec;
            outputArg2 = eta;
        end
    end
end

