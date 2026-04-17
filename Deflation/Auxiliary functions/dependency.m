function[check]=dependency(f, var)

% DEPENDENCY: checks anonymous function for variable dependency.
%
% check = DEPENDENCY(f,var) checks if the anonymous function f depends on
% the variable var (check = true) or not (check = false). The function f
% may be scalar valued or vector valued.
%
% Example:
% f = @(x,y) 2*x;
% check = dependency(f, 'x') is true.
% check = dependency(f, 'y') is false.

[~, remain]=strtok(func2str(f), ')');
check=contains(remain(2:end), var);
end