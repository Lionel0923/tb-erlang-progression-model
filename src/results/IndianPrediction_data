clear; clc; close all;

%% ==========================================================
% 1) Observed data (WHO 2010–2023)
%% ==========================================================
yearsObs = 2010:2023;

meanInc = [276 268 258 249 243 237 225 217 208 202 195 200 199 195]';
lowInc  = [114 138 141 141 171 181 191 193 185 180 172 173 170 164]';
uppInc  = [508 440 410 395 329 300 263 242 233 227 220 228 231 228]';

%% ==========================================================
% 2) Forecast data (2024–2033)
%% ==========================================================
yearsF = 2024:2033;

% ---- Simple baseline ----
simple_base_mean = [167.650 161.387 155.502 149.973 144.784 139.915 135.349 131.071 127.064 123.315]';
simple_base_lo   = [154.644 147.657 140.777 134.062 127.838 122.109 116.084 110.575 105.369 100.451]';
simple_base_hi   = [190.022 187.318 186.441 185.277 184.376 185.032 186.117 187.613 189.202 191.096]';

% ---- Erlang baseline ----
erlang_base_mean = [172.398 166.596 161.071 155.812 150.812 146.059 141.545 137.259 133.192 129.335]';
erlang_base_lo   = [159.194 152.694 146.435 140.196 134.422 128.590 123.188 118.408 113.848 109.467]';
erlang_base_hi   = [185.643 180.561 175.677 171.404 167.402 163.409 159.672 155.690 151.915 148.716]';

% ---- Simple ACF (2024–26) ----
simple_acf_mean = [149.491 134.832 124.955 129.706 128.599 124.271 120.260 116.545 113.084 109.858]';
simple_acf_lo   = [130.987 114.307 104.295 113.442 112.461 107.430 102.949  98.647  94.536  90.383]';
simple_acf_hi   = [170.114 156.359 146.224 148.211 152.736 154.554 156.386 158.445 160.823 163.563]';

% ---- Erlang ACF (2024–26) ----
erlang_acf_mean = [135.161 119.336 110.846 131.942 132.751 127.967 123.920 120.226 116.764 113.497]';
erlang_acf_lo   = [119.735  97.528  89.805 114.644 118.208 112.616 107.921 103.483  99.313  95.530]';
erlang_acf_hi   = [149.295 136.309 127.977 146.077 146.254 141.118 137.707 134.675 131.456 128.610]';

%% ==========================================================
% 3) Colors 
%% ==========================================================
simpleColor   = [0.00 0.35 0.90];   % Blue
erlangColor   = [0.90 0.10 0.10];   % Red
obsBandColor  = [0.78 0.86 0.62];   % Light green
obsMeanColor  = [0.20 0.50 0.10];   % Dark green

simpleBandColor = [0.55 0.70 1.00]; % Light blue band
erlangBandColor = [1.00 0.70 0.70]; % Light red band

%% ==========================================================
% 4) Plot
%% ==========================================================
figure('Position',[120 120 1100 620]); 
hold on; grid on; box on;

% --- ACF 95% UI bands ---
h_er_acf_ui = fill([yearsF fliplr(yearsF)], ...
                   [erlang_acf_lo' fliplr(erlang_acf_hi')], ...
                   erlangBandColor, ...
                   'EdgeColor','none','FaceAlpha',0.35, ...
                   'DisplayName','Erlang ACF 95% UI');

h_si_acf_ui = fill([yearsF fliplr(yearsF)], ...
                   [simple_acf_lo' fliplr(simple_acf_hi')], ...
                   simpleBandColor, ...
                   'EdgeColor','none','FaceAlpha',0.35, ...
                   'DisplayName','Simple ACF 95% UI');

% --- Observed 95% UI ---
h_obs_ui = fill([yearsObs fliplr(yearsObs)], ...
                [lowInc' fliplr(uppInc')], ...
                obsBandColor, ...
                'EdgeColor','none','FaceAlpha',0.45, ...
                'DisplayName','Observed 95% UI');

% --- Observed mean ---
h_obs_mean = plot(yearsObs, meanInc, 'o-', ...
                  'Color',obsMeanColor, ...
                  'LineWidth',3, ...
                  'MarkerFaceColor',obsMeanColor, ...
                  'DisplayName','Observed mean');

% --- Baseline lines ---
h_simple_base = plot(yearsF, simple_base_mean, '-', ...
                     'Color',simpleColor, ...
                     'LineWidth',3, ...
                     'DisplayName','Simple baseline');

h_erlang_base = plot(yearsF, erlang_base_mean, '-', ...
                     'Color',erlangColor, ...
                     'LineWidth',3, ...
                     'DisplayName','Erlang baseline');

% --- ACF dashed lines ---
h_simple_acf = plot(yearsF, simple_acf_mean, '--', ...
                    'Color',simpleColor, ...
                    'LineWidth',3, ...
                    'DisplayName','Simple (ACF 2024–26)');

h_erlang_acf = plot(yearsF, erlang_acf_mean, '--', ...
                    'Color',erlangColor, ...
                    'LineWidth',3, ...
                    'DisplayName','Erlang (ACF 2024–26)');

% ===== CONNECTORS: bridge 2023 -> 2024 (hidden from legend) =====
x23 = 2023;
x24 = 2024;
y23 = meanInc(end);  % observed mean in 2023 as anchor

% Baseline connectors (solid)
plot([x23 x24], [y23 simple_base_mean(1)], '-',  'Color', simpleColor, 'LineWidth',2, 'HandleVisibility','off');
plot([x23 x24], [y23 erlang_base_mean(1)], '-',  'Color', erlangColor, 'LineWidth',2, 'HandleVisibility','off');

% ACF connectors (dashed)
plot([x23 x24], [y23 simple_acf_mean(1)], '--', 'Color', simpleColor, 'LineWidth',2, 'HandleVisibility','off');
plot([x23 x24], [y23 erlang_acf_mean(1)], '--', 'Color', erlangColor, 'LineWidth',2, 'HandleVisibility','off');

%% ==========================================================
% 5) Axis & Legend
%% ==========================================================
xlabel('Year');
ylabel('Incidence per 100 000');
title('TB annual incidence in India');

xlim([2010 2033]);

allY = [lowInc; uppInc; ...
        simple_base_lo; simple_base_hi; ...
        erlang_base_lo; erlang_base_hi; ...
        simple_acf_lo; simple_acf_hi; ...
        erlang_acf_lo; erlang_acf_hi];

ylim([0 1.1*max(allY)]);

legend([h_er_acf_ui, h_si_acf_ui, h_obs_ui, h_obs_mean, ...
        h_simple_base, h_erlang_base, ...
        h_simple_acf, h_erlang_acf], ...
        'Location','northeast');

set(gca,'FontSize',12,'LineWidth',1);
