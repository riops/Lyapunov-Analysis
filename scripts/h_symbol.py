import json
from sympy import sqrt, Rational
from sympy.physics.wigner import clebsch_gordan, wigner_6j


def HSymbol(l1, m1, l2, m2, l3, m3, l4, m4, N):
    j = Rational(N - 1, 2)
    total = Rational(0)
    for l in range(N):
        p12 = 1 - (-1) ** (l1 + l2 + l)
        p34 = 1 - (-1) ** (l3 + l4 + l)
        six1 = wigner_6j(l1, l2, l, j, j, j)
        six2 = wigner_6j(l3, l4, l, j, j, j)
        for m in range(-l, l + 1):
            cg12 = clebsch_gordan(l1, l2, l, m1, m2, m)
            cg34 = clebsch_gordan(l3, l, l4, m3, m, m4)
            prefac = sqrt((2*l1 + 1)*(2*l2 + 1)*(2*l3 + 1)*(2*l + 1))
            term = prefac * p12 * p34 * six1 * six2 * cg12 * cg34
            total += term
    return total


def HContraction12(l1, m1, l2, m2, N):
    total = Rational(0)
    for l in range(N):
        for m in range(-l, l + 1):
            total += HSymbol(l, m, l, m, l2, m2, l1, m1, N)
    return total


def HContraction13(l1, m1, l2, m2, N):
    total = Rational(0)
    for l in range(N):
        for m in range(-l, l + 1):
            total += HSymbol(l, m, l2, m2, l, m, l1, m1, N)
    return total


def HContraction23(l1, m1, l2, m2, N):
    total = Rational(0)
    for l in range(N):
        for m in range(-l, l + 1):
            total += HSymbol(l2, m2, l, m, l, m, l1, m1, N)
    return total


if __name__ == '__main__':
    N_max = 4

    # 1) Grouped contractions into JSON by N
    # contraction_tasks = [
    #     (HContraction12, 'HContraction12_nonzero.json'),
    #     (HContraction13, 'HContraction13_nonzero.json'),
    #     (HContraction23, 'HContraction23_nonzero.json'),
    # ]
    #
    # for func, fname in contraction_tasks:
    #     grouped = {}
    #     for N in range(1, N_max + 1):
    #         entries = []
    #         print(f"Computing {fname} for N={N}...")
    #         for l1 in range(N):
    #             for m1 in range(-l1, l1 + 1):
    #                 for l2 in range(N):
    #                     for m2 in range(-l2, l2 + 1):
    #                         val = func(l1, m1, l2, m2, N)
    #                         if val != 0:
    #                             entries.append({
    #                                 'l1': l1,
    #                                 'm1': m1,
    #                                 'l2': l2,
    #                                 'm2': m2,
    #                                 'value': float(val.evalf())
    #                             })
    #         if entries:
    #             grouped[N] = entries
    #     with open(f'data/json/{fname}', 'w') as f:
    #         json.dump(grouped, f, indent=2)

    # 2) Combined sums grouped by N
    # combined_specs = [
    #     (lambda l1,m1,l2,m2,N: HContraction13(l1,m1,l2,m2,N) + HContraction23(l1,m1,l2,m2,N), 'HContraction13+HContraction23.json'),
    #     (lambda l1,m1,l2,m2,N: 2*HContraction13(l1,m1,l2,m2,N) + HContraction23(l1,m1,l2,m2,N), '2xHContraction13+HContraction23.json')
    # ]
    # for weight, fname in combined_specs:
    #     grouped = {}
    #     for N in range(1, N_max + 1):
    #         entries = []
    #         print(f"Computing {fname} for N={N}...")
    #         for l1 in range(N):
    #             for m1 in range(-l1, l1 + 1):
    #                 for l2 in range(N):
    #                     for m2 in range(-l2, l2 + 1):
    #                         val = weight(l1, m1, l2, m2, N)
    #                         if val != 0:
    #                             entries.append({
    #                                 'l1': l1,
    #                                 'm1': m1,
    #                                 'l2': l2,
    #                                 'm2': m2,
    #                                 'value': float(val.evalf())
    #                             })
    #         if entries:
    #             grouped[N] = entries
    #     with open(f'data/json/{fname}', 'w') as f:
    #         json.dump(grouped, f, indent=2)

    # 3) Nonzero HSymbol entries grouped by N
    grouped_symbol = {}
    fname = 'HSymbol_nonzero.json'
    for N in range(1, N_max + 1):
        entries = []
        print(f"Computing {fname} for N={N}...")
        for l1 in range(N):
            for m1 in range(-l1, l1 + 1):
                for l2 in range(N):
                    for m2 in range(-l2, l2 + 1):
                        for l3 in range(N):
                            for m3 in range(-l3, l3 + 1):
                                for l4 in range(N):
                                    for m4 in range(-l4, l4 + 1):
                                        val = HSymbol(l1, m1, l2, m2, l3, m3, l4, m4, N)
                                        if val != 0:
                                            entries.append({
                                                'l1': l1, 'm1': m1,
                                                'l2': l2, 'm2': m2,
                                                'l3': l3, 'm3': m3,
                                                'l4': l4, 'm4': m4,
                                                'value': float(val.evalf())
                                            })
        if entries:
            grouped_symbol[N] = entries
    with open(f'data/json/{fname}', 'w') as f:
        json.dump(grouped_symbol, f, indent=2)

    print("All grouped JSON files with floats have been written.")
