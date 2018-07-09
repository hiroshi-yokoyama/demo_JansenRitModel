function sigm = sigmoidal_function (e0, r, v0, v)
    sigm = 2*e0./(1 + exp(r .* (v0 - v) ) );
end