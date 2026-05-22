function Xay = enforce_age_year(X, nAge, nYear)
    if isequal(size(X), [nAge, nYear])
        Xay = X;
    elseif isequal(size(X), [nYear, nAge])
        Xay = X.';   
    else
        error('X must be %dx%d (age x year) or %dx%d (year x age). Got %dx%d.', ...
              nAge, nYear, nYear, nAge, size(X,1), size(X,2));
    end
end