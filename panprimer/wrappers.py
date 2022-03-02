import primer3


def design_primers(sequence, singleton=True, npairs=1, minlen=400, maxlen=4000):
    '''
    We don't have to explicitly exclude a region in the case of non-singletons
    bc/ we fill the space btw/ them w/ Ns and no Ns are allowed in primers.
    
    SEQUENCE_EXCLUDED_REGION would be the corresponding primer3 param.

    npairs = 1 .. best only

    The result of designPrimers() is ordered best to worse
    '''
    include = [0, len(sequence)]  # SEQUENCE_INCLUDED_REGION

    design = primer3.bindings.designPrimers(
        {
            'SEQUENCE_TEMPLATE': sequence,
            'SEQUENCE_INCLUDED_REGION': include,
        },
        {
            'PRIMER_NUM_RETURN': npairs,
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 24,
            'PRIMER_OPT_TM': 55.0,
            'PRIMER_MIN_TM': 50.0,
            'PRIMER_MAX_TM': 60.0,
            'PRIMER_MIN_GC': 40.0,
            'PRIMER_MAX_GC': 60.0,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_PRODUCT_SIZE_RANGE': [[minlen, maxlen]],
        })
    
    primers = []
    for i in range(npairs):
        if singleton:
            primers.append([
                design[f'PRIMER_LEFT_{i}_SEQUENCE'],
                design[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                round(design[f'PRIMER_LEFT_{i}_TM'], 2),
                round(design[f'PRIMER_RIGHT_{i}_TM'], 2),
                design[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
                round(design[f'PRIMER_PAIR_{i}_PENALTY'], 4),
                ])
        else:
            primers.append([
                design[f'PRIMER_LEFT_{i}_SEQUENCE'],
                design[f'PRIMER_RIGHT_{i}_SEQUENCE'],
                round(design[f'PRIMER_LEFT_{i}_TM'], 2),
                round(design[f'PRIMER_RIGHT_{i}_TM'], 2),
                'NA',
                round(design[f'PRIMER_PAIR_{i}_PENALTY'], 4),
                ])

    return primers