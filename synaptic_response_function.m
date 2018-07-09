function h = synaptic_response_function(A, a, t)
    h = zeros(size(t));
    h(t>=0) = A.*a.*t(t>=0).*exp(-a.*t(t>=0));
end