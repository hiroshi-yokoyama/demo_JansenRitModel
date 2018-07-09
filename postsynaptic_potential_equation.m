function f_out = postsynaptic_potential_equation(z, y, A, a, x)
    dy = z;
    dz = A.*a.*x - 2.*a.*z - a^2.*y; 
    f_out = [dy; dz];
end
