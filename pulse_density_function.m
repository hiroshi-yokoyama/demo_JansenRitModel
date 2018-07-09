function v_in = pulse_density_function (k)
    n = 7;
    w = 5;
    q = 0.0001;
    
    v_in = q.*(k./w).^n.*exp(-k./w);
end