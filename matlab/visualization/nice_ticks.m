function ticks = nice_ticks(maxVal)
    % Choose nice step size
    approxStep = maxVal / 10;
    pow10 = 10^floor(log10(approxStep));
    niceSteps = [1, 2, 5, 10];
    idx = find(niceSteps * pow10 >= approxStep, 1, 'first');
    if isempty(idx)
        % For very large step sizes, go up to the next power of 10
        step = 10 * pow10;
    else
        step = niceSteps(idx) * pow10;
    end

    % Create tick vector from 0 up to maxVal, step by 'step'
    ticks = 0:step:maxVal;
    % Make sure maxVal itself is included if it’s not exactly on a tick
    if ticks(end) < maxVal
        ticks(end+1) = ticks(end) + step;
    end
end
