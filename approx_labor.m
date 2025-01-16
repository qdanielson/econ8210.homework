function [l] = approx_labor(elements, weights)
    l = elements' * weights;
end