function updateProgress(total, reset)
% UPDATEPROGRESS displays calculation progress for a given cycle.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Displays information about the calculation progress for a given cycle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Miroslav HALLO
% Revision 2026/04: First version
%
% Copyright (C) 2026 Kyoto University
%
% This program is published under the GNU General Public License (GNU GPL).
%
% This program is free software: you can modify it and/or redistribute it
% or any derivative version under the terms of the GNU General Public
% License as published by the Free Software Foundation, either version 3
% of the License, or (at your option) any later version.
%
% This code is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY. We would like to kindly ask you to acknowledge the authors
% and don't remove their names from the code.
%
% You should have received copy of the GNU General Public License along
% with this program. If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%        total = total number of expected cycle executions
%        reset = true/false (not mandatory argument)
%
% OUTPUT:
%        None
%
% USAGE:
%        % Inicializace the progress function before the cycle
%        updateProgress([], true);
%        
%        for i = 1:N
%            % Display progress
%            updateProgress(N);
%            % ...
%            % It is NOT ALLOWED to display anything here (disp or fprintf)
%        end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Persistent variables
persistent count startTime;

% Reset progress function (clean persistent variables)
if nargin > 1 && reset == true
    count = [];
    startTime = [];
    return
end

% Not valid expected number of cycle executions
if total <= 0
    return
end

% First call
if isempty(count)
    count = 0;
    startTime = tic;
    fprintf('Processing:   0%% (ETA: --------)');
end

try
    % Compute percentage done
    count = count + 1;
    pct = floor((count/total) * 100);

    % Compute remaining time ETA
    if count > 3
        elapsedTime = toc(startTime);
        totalTime = (elapsedTime / (count-1)) * total;
        remainingTime = max(0, totalTime - elapsedTime);

        hours = floor(remainingTime / 3600);
        if hours < 100
            % Time format 'HH:MM:SS'
            mins = floor((remainingTime-(hours*3600)) / 60);
            secs = floor(mod((remainingTime-(hours*3600)), 60));
            etaStr = sprintf('%02d:%02d:%02d', hours, mins, secs);
        else
            % Time format 'DDD days'
            etaStr = sprintf('%03d days', round(hours/24));
        end
    else
        etaStr = '--------';
    end

    % Display progress
    fprintf(repmat('\b', 1, 20)); % Delete last chars
    fprintf('%3d%% (ETA: %s)', pct, etaStr); % Display new

    % End progress and clean persistent variables
    if count >= total
        % Total elapsed time
        elapsedTime = toc(startTime);
        hours = floor(elapsedTime / 3600);
        if hours < 100
            % Time format 'HH:MM:SS'
            mins = floor((elapsedTime-(hours*3600)) / 60);
            secs = floor(mod((elapsedTime-(hours*3600)), 60));
            elaStr = sprintf('%02d:%02d:%02d', hours, mins, secs);
        else
            % Time format 'DDD days'
            elaStr = sprintf('%03d days', round(hours/24));
        end
        fprintf(repmat('\b', 1, 15));  % Delete last chars
        fprintf('Done (Elapsed time: %s).\n', elaStr); % Display new
        count = [];
        startTime = [];
    end

    drawnow('limitrate');

catch
    % In the case of unexpected error
    fprintf(repmat('\b', 1, 20));
    fprintf('Wait|No ETA estimate');
end

end

