
function [ ] = checkPrompt( value, condition )
% Checks if the prompted value satisfies the condition

if isempty(value)
    error('ERROR: please enter a value.')
end

switch condition
    case 'method'
        if value < 0 || value > 3 || mod(value, 1) ~= 0
            error('ERROR: invalid method index.')
        end
    case 'nonnegInteger'
        if value < 0 || mod(value, 1) ~= 0
            error('ERROR: please enter a positive integer.')
        end
    case '>=1Integer'
        if value < 1 || mod(value, 1) ~= 0
            error('ERROR: please enter an integer >= 1.')
        end
end

end

