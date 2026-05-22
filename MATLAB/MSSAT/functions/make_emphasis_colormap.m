function [cmapAge, t_of_age] = make_emphasis_colormap(ages, emphBand, baseMap, tEmph)
    % ages: 0..100 (or any sorted vector)
    % emphBand: [aLo aHi], ages in this band get mapped to tEmph interval in base colormap
    % baseMap: Mx3
    % tEmph: [tLo tHi] subset of [0,1] where colors are more salient

    aMin = min(ages); aMax = max(ages);
    aLo  = emphBand(1); aHi = emphBand(2);

    tLo = tEmph(1); tHi = tEmph(2);
    t_of_age = zeros(size(ages));

    for i = 1:numel(ages)
        a = ages(i);
        if a <= aLo
            % [aMin, aLo] -> [0, tLo]
            if aLo == aMin
                t = 0;
            else
                t = (a - aMin) / (aLo - aMin) * tLo;
            end
        elseif a <= aHi
            % [aLo, aHi] -> [tLo, tHi]  
            if aHi == aLo
                t = (tLo + tHi)/2;
            else
                t = tLo + (a - aLo) / (aHi - aLo) * (tHi - tLo);
            end
        else
            % [aHi, aMax] -> [tHi, 1]
            if aMax == aHi
                t = 1;
            else
                t = tHi + (a - aHi) / (aMax - aHi) * (1 - tHi);
            end
        end
        t_of_age(i) = min(max(t,0),1);
    end

    % sample baseMap at these t values
    M = size(baseMap,1);
    idx = 1 + round(t_of_age * (M-1));
    idx = min(max(idx,1),M);
    cmapAge = baseMap(idx, :);
end