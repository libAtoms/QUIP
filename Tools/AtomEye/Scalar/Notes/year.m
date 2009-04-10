%%%%%%%%%%
% year.m %
%%%%%%%%%%

% e = 1/2, 1 period %

load year_data.out;
methods = 13;
n = size(year_data,1);
i = 1;

Ruth83 = year_data(i:methods:n,:); i=i+1;
Schlier98_6a = year_data(i:methods:n,:); i=i+1;
Tsitouras99 = year_data(i:methods:n,:); i=i+1;

Calvo93 = year_data(i:methods:n,:); i=i+1;
Schlier00_6b = year_data(i:methods:n,:); i=i+1;
Schlier00_8b = year_data(i:methods:n,:); i=i+1;
Schlier00_8c = year_data(i:methods:n,:); i=i+1;

rk4 = year_data(i:methods:n,:); i=i+1;
gear4 = year_data(i:methods:n,:); i=i+1;
gear5 = year_data(i:methods:n,:); i=i+1;
gear6 = year_data(i:methods:n,:); i=i+1;
gear7 = year_data(i:methods:n,:); i=i+1;
gear8 = year_data(i:methods:n,:); i=i+1;

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

title ([ 'Integration of 1 period of ' ...
         'Kepler orbital with eccentricity 0.5' ]);
set (gca, 'XTick', [50 100:100:1000] );
axis ([50 1000 1e-11 5]);
print -depsc year_a.eps;

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

title ([ 'Integration of 1 period of ' ...
         'Kepler orbital with eccentricity 0.5' ]);
set (gca, 'XTick', [50 100:100:1000] );
axis ([50 1000 1e-16 1]);

print -depsc year_b.eps;
