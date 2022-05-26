function b = isValidPositiveScalar(x)
if (isnumeric(x) && isscalar(x) && (x > 0))
    b = true;
else
    b = false;
end
end
