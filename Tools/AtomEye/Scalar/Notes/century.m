%%%%%%%%%%%%%
% century.m %
%%%%%%%%%%%%%

% e = 1/2, 100 periods %

load century_data.out;
methods = 13;
n = size(century_data,1);
i = 1;

Ruth83 = century_data(i:methods:n,:); i=i+1;
Schlier98_6a = century_data(i:methods:n,:); i=i+1;
Tsitouras99 = century_data(i:methods:n,:); i=i+1;

Calvo93 = century_data(i:methods:n,:); i=i+1;
Schlier00_6b = century_data(i:methods:n,:); i=i+1;
Schlier00_8b = century_data(i:methods:n,:); i=i+1;
Schlier00_8c = century_data(i:methods:n,:); i=i+1;

rk4 = century_data(i:methods:n,:); i=i+1;
gear4 = century_data(i:methods:n,:); i=i+1;
gear5 = century_data(i:methods:n,:); i=i+1;
gear6 = century_data(i:methods:n,:); i=i+1;
gear7 = century_data(i:methods:n,:); i=i+1;
gear8 = century_data(i:methods:n,:); i=i+1;

clf;
loglog (Ruth83(:,1), Ruth83(:,2), 'go');
hold on;
loglog (Schlier98_6a(:,1), Schlier98_6a(:,2), 'gs');
loglog (Tsitouras99(:,1), Tsitouras99(:,2), 'gx');

loglog (Calvo93(:,1), Calvo93(:,2), 'md');
loglog (Schlier00_6b(:,1), Schlier00_6b(:,2), 'r*');
loglog (Schlier00_8c(:,1), Schlier00_8c(:,2), 'k^');

loglog (rk4(:,1), rk4(:,2), '+');
loglog (gear4(:,1), gear4(:,2), '>');
loglog (gear5(:,1), gear5(:,2), 'c<');
loglog (gear6(:,1), gear6(:,2), 'p');
loglog (gear7(:,1), gear7(:,2), 'ch');
loglog (gear8(:,1), gear8(:,2), 'v');

legend ('Ruth83', 'Schlier98\_6a', 'Tsitouras99', ...
        'Calvo93', 'Schlier00\_6b', 'Schlier00\_8c', ...
        '4th Runge-Kutta', '4th Gear', '5th Gear', ...
        '6th Gear', '7th Gear', '8th Gear', 3);

loglog (Schlier00_6b(:,1), Schlier00_6b(:,2), 'r');

xlabel ('number of force evaluations per period');
sylabel ('|| final (\bold p\normal,\bold q\normal) error ||_2');

title ([ 'Integration of 100 periods of ' ...
         'Kepler orbitals with eccentricity 0.5' ]);
set (gca, 'XTick', [100 150 200:100:1000] );
axis ([100 1000 1e-6 5]);
print -depsc century_a.eps;

input('Press a key... ');

clf;
loglog (Ruth83(:,1), abs(Ruth83(:,3)), 'go');
hold on;
loglog (Schlier98_6a(:,1), abs(Schlier98_6a(:,3)), 'gs');
loglog (Tsitouras99(:,1), abs(Tsitouras99(:,3)), 'gx');

loglog (Calvo93(:,1), abs(Calvo93(:,3)), 'md');
loglog (Schlier00_6b(:,1), abs(Schlier00_6b(:,3)), 'r*');
loglog (Schlier00_8c(:,1), abs(Schlier00_8c(:,3)), 'k^');

loglog (rk4(:,1), abs(rk4(:,3)), '+');
loglog (gear4(:,1), abs(gear4(:,3)), '>');
loglog (gear5(:,1), abs(gear5(:,3)), 'c<');
loglog (gear6(:,1), abs(gear6(:,2)), 'p');
loglog (gear7(:,1), abs(gear7(:,2)), 'ch');
loglog (gear8(:,1), abs(gear8(:,2)), 'v');

legend ('Ruth83', 'Schlier98\_6a', 'Tsitouras99', ...
        'Calvo93', 'Schlier00\_6b', 'Schlier00\_8c', ...
        '4th Runge-Kutta', '4th Gear', '5th Gear', ...
        '6th Gear', '7th Gear', '8th Gear', 3);

loglog (Schlier00_6b(:,1), abs(Schlier00_6b(:,3)), 'r');

xlabel ('number of force evaluations per period');
ylabel ('| final energy error |');

title ([ 'Integration of 100 periods of ' ...
         'Kepler orbitals with eccentricity 0.5' ]);
set (gca, 'XTick', [100 150 200:100:1000] );
axis ([100 1000 1e-15 1]);

print -depsc century_b.eps;
