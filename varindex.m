function ind = varindex(terms, variables)
    if ~iscell(terms)
        terms = {terms};
    end
    m = length(terms);
    n = length(variables);
    l = cellfun(@length, variables);
    assert(all(l == l(1)));
    l = l(1);
    ind = zeros(1, m);
    f = [1 cumprod(ones(1, n-1)*2)];
    for i = 1:m
        t = terms{i};
        sub = zeros(1, n);
        while ~isempty(t)
            [~, j] = ismember(t(1:l), variables);
            assert(length(j) == 1, 'variable %s not found', t(1:l));
            assert(sub(j) == 0, 'variable %s is present more than once', t(1:l));
            sub(j) = 1;
            t = t(l+1:end);
        end
        ind(i) = sum(sub.*f) + 1;
    end
end
