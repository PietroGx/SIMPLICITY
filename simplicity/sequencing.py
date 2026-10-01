# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# Sequencing datasets, built AFTER the simulation rather than during it.
# ----------------------------------------------------------------------------
# The simulation used to draw a sequencing sample inline, on every diagnosis,
# against `sequencing_rate`. That mixed a surveillance choice into the run: a
# ratio could only be changed by simulating again, and at the production rate
# of 0.05 one seed yielded ~7 sequences out of ~1,500 infections.
#
# Now a run records `t_diagnosis` on every diagnosed individual and writes the
# COMPLETE record -- every diagnosed host, every distinct genome it carried.
# That is what rate 1.0 produced before, which is what both calibration stages
# already asked for. Any lower ratio is a subsample of that record, drawn here
# from saved output with `reconstruct_sequencing_data`, as many times and at as
# many ratios as you like.
#
# A host is sequenced or not, and then ALL of its distinct genomes are read:
# the draw is per individual, never per row, matching the mechanism this
# replaces.
# ============================================================================

import ast
import re

import numpy as np

ROW_FIELDS = ('individual_index', 'sequencing_time', 'lineage_name',
              'sequence', 'sequence_lenght', 'individual_type',
              'infection_duration', 'intra-host_lineage_index')


def _distinct_genomes(lineage_names, genome_of):
    '''(lineage_name, genome) per DISTINCT genome, in carriage order.

    IH_lineages is a multiset: add_lineage duplicates an existing lineage into
    a new slot and that copy stays byte-identical until it mutates, so emitting
    every slot would record the same genome several times and weight the host
    by its lineage count instead of its actual diversity.
    '''
    seen = set()
    out = []
    for lineage_name in lineage_names:
        genome = genome_of(lineage_name)
        if genome is None:
            continue
        key = tuple(sorted(genome.items()))
        if key in seen:
            continue
        seen.add(key)
        out.append((lineage_name, genome))
    return out


def _row(ind_index, t_diagnosis, lineage_name, genome, ind_type, t_infection,
         ih_index):
    return {
        'individual_index': ind_index,
        'sequencing_time': t_diagnosis,
        'lineage_name': lineage_name,
        'sequence': genome,
        'sequence_lenght': len(genome),
        'individual_type': ind_type,
        'infection_duration': t_diagnosis - t_infection,
        'intra-host_lineage_index': ih_index,
    }


def complete_sequencing_rows(population):
    '''Every diagnosed host's distinct genomes, from the live population.

    Called at save time, so it sees the whole run. A host leaves infected_i on
    diagnosis and stops acquiring or mutating lineages there, so its
    IH_lineages at the end of the run are the ones it carried when diagnosed.

    Ordered by diagnosis time, so the file reads chronologically exactly as the
    inline draw's did.
    '''
    diagnosed = [(ind['t_diagnosis'], i, ind)
                 for i, ind in population.individuals.items()
                 if ind.get('t_diagnosis') is not None]
    diagnosed.sort(key=lambda t: (t[0], t[1]))

    rows = []
    for t_diagnosis, i, ind in diagnosed:
        pairs = _distinct_genomes(ind['IH_lineages'],
                                  population.get_lineage_genome)
        for ih_index, (lineage_name, genome) in enumerate(pairs):
            rows.append(_row(i, t_diagnosis, lineage_name, genome,
                             ind['type'], ind['t_infection'], ih_index))
    return rows


def subsample(rows, ratio, seed=None):
    '''Keep each diagnosed INDIVIDUAL with probability `ratio`.

    Per individual, not per row: a host is sequenced or it is not, and then all
    of its distinct genomes are read. Subsampling rows independently would
    break the host-level sampling the surveillance model represents.
    '''
    if not 0.0 <= ratio <= 1.0:
        raise ValueError(f'ratio must be in [0, 1], got {ratio!r}')
    if ratio == 1.0:
        return list(rows)
    if ratio == 0.0:
        return []

    rng = np.random.default_rng(seed)
    # one draw per individual, in first-appearance order, so the result is
    # reproducible from (rows, ratio, seed) alone
    keep = {}
    for row in rows:
        i = row['individual_index']
        if i not in keep:
            keep[i] = rng.uniform() < ratio
    return [row for row in rows if keep[row['individual_index']]]


_NP_SCALAR = re.compile(r'np\.(?:float64|float32|int64|int32|bool_)\(([^()]*)\)')


def _literal(value, default):
    '''ast.literal_eval tolerant of the numpy scalar reprs these CSVs carry.'''
    if not isinstance(value, str):
        return default if value is None else value
    try:
        return ast.literal_eval(_NP_SCALAR.sub(r'\1', value))
    except (ValueError, SyntaxError):
        return default


def rows_from_saved_output(individuals_data, phylogenetic_data):
    '''The complete record, rebuilt from a finished run's CSVs.

    individuals_data -- as read_individuals_data returns it, carrying
                        t_diagnosis, IH_lineages, type and t_infection
    phylogenetic_data -- as read_phylogenetic_data returns it, carrying
                        Lineage_name and Genome

    Use this to regenerate a sequencing dataset for a run that is already on
    disk; `complete_sequencing_rows` is the same thing computed in-process.
    '''
    if 't_diagnosis' not in individuals_data.columns:
        raise KeyError(
            "individuals_data has no 't_diagnosis' column: it predates the "
            "removal of the in-simulation sequencing draw, so a sequencing "
            "dataset cannot be reconstructed from it. Use its own "
            "sequencing_data.csv, or rerun.")

    genomes = {}
    for _, prow in phylogenetic_data.iterrows():
        genomes[prow['Lineage_name']] = _literal(prow['Genome'], {})

    diagnosed = []
    for i, irow in individuals_data.iterrows():
        t_diagnosis = irow['t_diagnosis']
        if t_diagnosis is None or (isinstance(t_diagnosis, float)
                                   and np.isnan(t_diagnosis)):
            continue
        diagnosed.append((float(t_diagnosis), i, irow))
    diagnosed.sort(key=lambda t: (t[0], t[1]))

    rows = []
    for t_diagnosis, i, irow in diagnosed:
        lineage_names = _literal(irow['IH_lineages'], [])
        pairs = _distinct_genomes(lineage_names, genomes.get)
        for ih_index, (lineage_name, genome) in enumerate(pairs):
            rows.append(_row(i, t_diagnosis, lineage_name, genome,
                             irow['type'], float(irow['t_infection']),
                             ih_index))
    return rows


def reconstruct_sequencing_data(individuals_data, phylogenetic_data,
                                sequencing_ratio=1.0, seed=None):
    '''A sequencing dataset at any ratio, from a run already on disk.

    This is the post-hoc replacement for the old `sequencing_rate` parameter:
    the ratio is now an argument you can vary without simulating again.

        import simplicity.output_manager as om
        import simplicity.sequencing as seq
        rows = seq.reconstruct_sequencing_data(
            om.read_individuals_data(ssod), om.read_phylogenetic_data(ssod),
            sequencing_ratio=0.05, seed=1)
    '''
    return subsample(rows_from_saved_output(individuals_data,
                                            phylogenetic_data),
                     sequencing_ratio, seed)
