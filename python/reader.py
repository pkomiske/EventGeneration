import time
import warnings

import numpy as np

allowed_particle_codes = frozenset(['C', 'R', 'O', 'D'])

#function to read in generated jets
def read(filepath, num=-1, max_files=None, verbose=1,
                  store_hadrons=True, store_partons=True, 
                  which_particles='all'):

    # start timing
    start = time.time()

    # check for properly formed options
    if which_particles is not None:
        which_particles = which_particles.upper()

    if which_particles == 'ALL':
        store_all = store_origs = store_decays = True
        store_C_only = store_R_only = False

    elif which_particles == 'NONE' or which_particles is None:
        store_all = store_origs = store_decays = store_C_only = store_R_only = False

    else:

        # check all letters understood
        for letter in which_particles:
            if letter not in allowed_particle_codes:
                raise ValueError("'{}' not in allowed particle codes: {}".format(letter, allowed_particle_codes))

        store_origs = ('O' in which_particles)
        store_decays = ('D' in which_particles)
        store_all = ('C' in which_particles and 'R' in which_particles)
        store_C_only = ((not store_all) and 'C' in which_particles)
        store_R_only = ((not store_all) and 'R' in which_particles)

    # derived booleans
    store_any = (store_all or store_C_only or store_R_only)
    store_hadron_particles = (store_hadrons and store_any)
    store_parton_particles = (store_partons and store_any)

    # arrays
    statuses, weights = [], []
    origs_h, origs_p = [], []
    decays_h, decays_p = [], []
    jets_h, jets_p = [], []
    hadrons, partons = [], []
    Cmasks_h, Cmasks_p = [], []
    Rmasks_h, Rmasks_p = [], []
    all_mask_bits = []
    sigmas, sigma_errs = [], []
    total_weights = []
    durations = []

    # filepath is a pattern and we utilize max_files
    if max_files is not None:
        filepaths = [filepath.format(i) for i in range(max_files)]

    # filepath is a single file
    elif isinstance(filepath, str):
        filepaths = [filepath]

    # filepath is some container of filepaths
    else:
        filepaths = filepath

    # iterate over files and count how many total units we have
    ntot = 0
    for file_i,filepath in enumerate(filepaths):
        with open(filepath, 'r') as f:

            # consume comments
            comments = []
            all_mask_bits.append({})
            for row in f:
                if not row.startswith('#'):
                    if verbose >= 2:
                        print(''.join(comments))

                    # check consistent mask bits
                    if file_i > 0 and all_mask_bits[0] != all_mask_bits[-1]:
                        m = 'Mask bits differ from rest in file {}: {}\nvs.\n{}'
                        warnings.warn(m.format(filepath, all_mask_bits[0], all_mask_bits[-1]))

                    # exit comments loop
                    break

                # store comments
                comments.append(row)

                # store mask bits
                if 'M_' in row:
                    name, val = row[1:].split(':')
                    all_mask_bits[-1][name.strip()] = int(val)

            # iterate over units in file
            fcount = 0
            for row in f:

                parts = row.split()
                if len(parts) == 0:

                    # check ending condition
                    if ntot == num:
                        break

                    # start a new jet with every newline
                    if store_hadrons:
                        if store_origs:
                            origs_h.append(orig_h)
                        if store_decays:
                            decays_h.append(decay_h)
                        jets_h.append(jet_h)
                        if store_any:
                            hadrons.append(np.asarray(hadrons_i, dtype=float))
                            if store_all:
                                Cmasks_h.append(np.asarray(Cmask_h, dtype=bool))
                                Rmasks_h.append(np.asarray(Rmask_h, dtype=bool))

                    if store_partons:
                        if store_origs:
                            origs_p.append(orig_p)
                        if store_decays:
                            decays_p.append(decay_p)
                        jets_p.append(jet_p)
                        if store_any:
                            partons.append(np.asarray(partons_i, dtype=float))
                            if store_all:
                                Cmasks_p.append(np.asarray(Cmask_p, dtype=bool))
                                Rmasks_p.append(np.asarray(Rmask_p, dtype=bool))

                    ntot += 1
                    fcount += 1
                    continue

                # get the first character of the key that starts the line
                key = parts[0]
                key0 = key[0]

                # hadron
                if key0 == 'H':
                    if store_hadron_particles:
                        if store_all:
                            hadrons_i.append(parts[1:])
                            Cmask_h.append('C' in key)
                            Rmask_h.append('R' in key)

                        elif store_C_only:
                            if 'C' in key:
                                hadrons_i.append(parts[1:])

                        elif store_R_only:
                            if 'R' in key:
                                hadrons_i.append(parts[1:])
                    continue

                # parton
                if key0 == 'P':
                    if store_parton_particles:
                        if store_all:
                            partons_i.append(parts[1:])
                            Cmask_p.append('C' in key)
                            Rmask_p.append('R' in key)

                        elif store_C_only:
                            if 'C' in key:
                                partons_i.append(parts[1:])

                        elif store_R_only:
                            if 'R' in key:
                                partons_i.append(parts[1:])
                    continue

                # decayed particle
                if key0 == 'D':
                    if store_decays:
                        decay.append(parts[1:])
                    continue

                # jet
                if key0 == 'J':
                    key1 = key[1]

                    # hadron jet
                    if key1 == 'H':
                        if store_hadrons:
                            jet_h = np.asarray(parts[1:], dtype=float)
                        continue

                    # parton jet
                    if key1 == 'P':
                        if store_partons:
                            jet_p = np.asarray(parts[1:], dtype=float)
                        continue

                    raise RuntimeError("Unknown line '" + key + "'")

                # original (hard-process) particle
                if key0 == 'O':
                    if store_origs:
                        if key == 'O':
                            orig_p = orig_h = np.asarray(parts[1:], dtype=float)
                            decay_h = decay_p = decay = []
                        elif 'H' in key:
                            orig_h = np.asarray(parts[1:], dtype=float)
                            decay_h = decay = []
                        elif 'P' in key:
                            orig_p = np.asarray(parts[1:], dtype=float)
                            decay_p = decay = []
                        else:
                            raise RuntimeError("Unknown key '" + key + "'")
                    continue

                # event status
                if key0 == 'S':
                    statuses.append(parts[1])
                    continue

                # event weight
                if key0 == 'W':
                    weights.append(parts[1])
                    continue

                # unit index, reset event containers
                if key0 == 'i':
                    if store_hadrons:
                        hadrons_i, Cmask_h, Rmask_h = [], [], []
                        jet_h = orig_h = decay_h = None

                    if store_partons:
                        partons_i, Cmask_p, Rmask_p = [], [], []
                        jet_p = orig_p = decay_p = None

                    continue

                # ending comments
                if key0 == '#':
                    key1 = parts[1] if len(parts) > 1 else ''
                    if key1.startswith('CrossSection'):
                        sigmas.append(parts[2])
                        sigma_errs.append(parts[4])

                    elif key1.startswith('TotalWeight'):
                        total_weights.append(parts[2])

                    elif key1.startswith('Duration'):
                        durations.append(parts[2][:-1])

                    continue

                # unknown key
                raise RuntimeError("Unknown key '" + key + "'")

        # print update
        if verbose >= 1:
            s = '\n' if verbose >= 2 else ''
            print('{:.3f}s - Read {} units from file {}{}'.format(time.time() - start, fcount, filepath, s))

        # break out of loop over files
        if ntot == num:
            break

    # dictionary to hold all the arrays
    d = {'mask_bits': all_mask_bits[0],
         'statuses': np.asarray(statuses, dtype=int), 
         'weights': np.asarray(weights, dtype=float),
         'sigmas': np.asarray(sigmas, dtype=float),
         'sigma_errs': np.asarray(sigma_errs, dtype=float),
         'total_weights': np.asarray(total_weights, dtype=float),
         'durations': np.asarray(durations, dtype=float)}

    if store_hadrons:
        d['origs_h'] = np.asarray(origs_h, dtype=object)
        d['decays_h'] = np.asarray(decays_h, dtype=object)
        d['jets_h'] = np.asarray(jets_h, dtype=object)
        d['hadrons'] = np.asarray(hadrons, dtype=object)
        if store_all:
            d['Cmasks_h'] = np.asarray(Cmasks_h, dtype=object)
            d['Rmasks_h'] = np.asarray(Rmasks_h, dtype=object)

    if store_partons:
        d['origs_p'] = np.asarray(origs_p, dtype=object)
        d['decays_p'] = np.asarray(decays_p, dtype=object)
        d['jets_p'] = np.asarray(jets_p, dtype=object)
        d['partons'] = np.asarray(partons, dtype=object)
        if store_all:
            d['Cmasks_p'] = np.asarray(Cmasks_p, dtype=object)
            d['Rmasks_p'] = np.asarray(Rmasks_p, dtype=object)

    return d
