function O = compileE1E2(settings,E1,E2)
    if nargin < 1 || isempty(settings)
        settings = sdpsettings('solver', 'linprog');
    end

    sdpvar E1 E2

    vO = {'A' 'B' 'C' 'D'};
    vI = {'A1' 'B1' 'C1' 'D1' 'A2' 'B2'};

    constraints = {'A1B1C1D1' 'ABCD' % new set
                   'A1B1C1' 'ABC'
                   'A1B1D1' 'ABD'
                   'A1C1D1' 'ACD'
                   'B1C1D1' 'BCD'
                   'A1B1' 'AB'
                   'B1C1' 'BC'
                   'C1D1' 'CD'
                   'A1D1' 'AD'
                   'A1C1' 'A*C'
                   'B1D1' 'B*D'
                   'A1' 'A'
                   'B1' 'B'
                   'C1' 'C'
                   'D1' 'D'
                   '' ''
                   'A2B2C1D1' 'ABCD' % new set
                   'A2B2C1' 'ABC'
                   'A2B2D1' 'ABD'
                   'A2C1D1' 'ACD'
                   'B2C1D1' 'BCD'
                   'A2B2' 'AB'
                   'B2C1' 'BC'
                   'A2D1' 'AD'
                   'A2C1' 'A*C'
                   'B2D1' 'B*D'
                   'A1B2D1' 'AD*B' % new set
                   'A1B2' 'A*B'
                   'A2B1C1' 'A*BC' % new set
                   'A2B1' 'A*B'
                   'A2' 'A'
                   'B2' 'B'};

    symConstraints = {{'A1A2B1' 'A1A2B2'}
                      {'A1B1B2' 'A2B1B2'}};

    nO = length(vO);
    nI = length(vI);

    CO = sdpvar(2^nO, 1);
    CI = sdpvar(2^nI, 1);
    slack = sdpvar;
    CO(varindex('', vO)) = 1;
    CO(varindex({'A' 'B' 'C' 'D'}, vO)) = E1;
    CO(varindex({'AB' 'BC' 'CD' 'AD'}, vO)) = E2;
    CO(varindex({'AC' 'BD'}, vO)) = E1*E1;

    cnames = constraints(:,1);
    snames = {};
    for i = 1:length(symConstraints)
        snames = horzcat(snames, symConstraints{i});
    end
    assert(length(unique(cnames)) == length(cnames), 'Duplicate constraints');
    assert(isempty(intersect(cnames, snames)), 'constraints and symConstraints should not overlap');
    assert(length(unique(snames)) == length(snames), 'Duplcate symConstraints');

    F = [];
    for i = 1:size(constraints, 1)
        lhs = varindex(constraints{i, 1}, vI);
        rhs = varindex(strsplit(constraints{i, 2}, '*'), vO);
        lhs = CI(lhs);
        rhs = prod(CO(rhs));
        F = [F; lhs == rhs];
    end
    H = [1 1
         1 -1];
    M = 1;
    for i = 1:nI
        M = kron(M, H);
    end
    F = [F; M * CI >= slack];
    O = optimizer(F, -slack, settings, {E1, E2}, slack);
end
