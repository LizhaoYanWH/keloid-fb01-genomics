"""Figures 1 and 2 reconstructed only from SHA-bound completed source tables.

No code runs on import. No model, inference, UMAP fit, scoring, sampling or FDR.
Expected ctx: data (source-id -> DataFrame), save(fig,stem,caption,sourcekeys),
panel(ax,letter,title). Root's common context owns backend/export/provenance.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

BLUE = '#0072B2'
ORANGE = '#D55E00'
GREEN = '#009E73'
GRAY = '#666666'
COHORTS = ['GSE241132_D0', 'GSE220300', 'GSE335482']
COHORT_COLORS = {'GSE241132_D0': GRAY, 'GSE220300': BLUE,
                 'GSE335482': ORANGE, 'GSE163973': GREEN}
LINEAGE_ORDER = ['fibroblast', 'keratinocyte', 'endothelial',
                 'lymphatic_endothelial', 'mural', 'myeloid', 'lymphoid',
                 't_nk', 'mast', 'melanocyte', 'schwann',
                 'glandular_epithelial', 'erythroid', 'unknown']
LINEAGE_LABELS = {'fibroblast': 'Fibroblast', 'keratinocyte': 'Keratinocyte',
    'endothelial': 'Endothelial', 'lymphatic_endothelial': 'Lymphatic endo.',
    'mural': 'Mural', 'myeloid': 'Myeloid', 'lymphoid': 'Lymphoid',
    't_nk': 'T/NK', 'mast': 'Mast', 'melanocyte': 'Melanocyte',
    'schwann': 'Schwann', 'glandular_epithelial': 'Gland epith.',
    'erythroid': 'Erythroid', 'unknown': 'Unknown'}
LINEAGE_PALETTE = ['#0072B2', '#E69F00', '#009E73', '#56B4E9', '#CC79A7',
    '#D55E00', '#7851A9', '#AA4499', '#A6761D', '#999933', '#44AA99',
    '#88CCEE', '#882255', '#B5B5B5']
STATE_PALETTE = ['#A6CEE3', '#0072B2', '#B2DF8A', '#FDBF6F', '#CAB2D6',
                 '#FB9A99', '#80B1D3', '#B3B3B3', '#B3DE69']


def _frame(ctx, key, columns, rows=None):
    if key not in ctx.data:
        raise KeyError('Required completed source: ' + key)
    frame = ctx.data[key]
    absent = set(columns).difference(frame.columns)
    if absent:
        raise ValueError(key + ' missing columns: ' + str(sorted(absent)))
    if rows is not None and len(frame) != rows:
        raise ValueError(key + ' wrong complete-family row count')
    return frame.copy()


def _numeric(values):
    """Decode original numerical columns; only explicit missing values are NA."""
    result = pd.to_numeric(values, errors='coerce')
    originals = pd.Series(values)
    missing = originals.isna() | originals.astype(str).str.strip().str.lower().isin(
        ['', 'na', 'nan', 'none'])
    if bool((pd.isna(result) & ~missing.to_numpy()).any()):
        raise ValueError('Unexpected non-numeric source value; do not silently erase')
    return np.asarray(result, dtype=float)


def _axis(fig, rect):
    ax = fig.add_axes(rect)
    ax.tick_params(labelsize=5.5, width=.5, length=2)
    for spine in ax.spines.values():
        spine.set_linewidth(.5)
    return ax


def _umap(ax, frame, category, colors, point_size=.38):
    if not frame.cell_id.is_unique:
        raise ValueError('Original cell IDs must remain unique')
    if not np.isfinite(frame[['UMAP1', 'UMAP2']].to_numpy(float)).all():
        raise ValueError('Original coordinates must be complete and finite')
    categories = frame[category].astype(str)
    if set(categories).difference(colors):
        raise ValueError('Unregistered original coordinate category')
    for label in colors:
        subset = frame.loc[categories.eq(label)]
        if len(subset):
            ax.scatter(subset.UMAP1, subset.UMAP2, s=point_size,
                       c=colors[label], linewidths=0, alpha=.82,
                       rasterized=True, marker='o')
    ax.set_xlabel('UMAP1', fontsize=6.5, labelpad=1)
    ax.set_ylabel('UMAP2', fontsize=6.5, labelpad=1)
    ax.xaxis.set_major_locator(MaxNLocator(3))
    ax.yaxis.set_major_locator(MaxNLocator(3))
    ax.set_aspect('equal', adjustable='datalim')
    ax.margins(.04)


def _dot_matrix(ax, data, family, row_order, columns, norm, cmap):
    selected = data.loc[data.group_family.eq(family)].copy()
    # Case normalization is for display only. Lymphoid and T_NK remain separate.
    selected['display_row'] = (selected.group.astype(str).str.lower()
                               if family == 'lineage' else selected.group.astype(str))
    if selected[['display_row', 'cohort']].duplicated().any():
        raise ValueError('Display normalization must not collapse source tests')
    seen = set()
    for row in selected.itertuples():
        iy = row_order.index(row.display_row)
        ix = columns.index(row.cohort)
        q = float(row.FDR_BH_all47_fixed_MC)
        if not np.isfinite(q) or not 0 < q <= 1:
            raise ValueError('Invalid existing q value')
        area = 6.0 * (1.0 + np.log10(max(float(row.n_cells), 1.0)))
        ax.scatter(ix, iy, s=area, c=[cmap(norm(-np.log10(q)))],
                   edgecolors='black' if q < .05 else '#888888',
                   linewidths=.8 if q < .05 else .35)
        seen.add((ix, iy))
    # A gray cross is a group not present in that original cohort taxonomy.
    for iy in range(len(row_order)):
        for ix in range(len(columns)):
            if (ix, iy) not in seen:
                ax.scatter(ix, iy, marker='x', s=8, c='#B6B6B6', linewidths=.5)
    ax.set_xlim(-.55, len(columns)-.45)
    ax.set_ylim(len(row_order)-.5, -.5)
    ax.set_xticks(range(len(columns)),
                  [dict(zip(COHORTS, ['D0', '220', '335']))[x] for x in columns])
    ax.set_yticks(range(len(row_order)),
                  [LINEAGE_LABELS.get(x, x) for x in row_order])
    ax.tick_params(axis='both', length=0, labelsize=5.5)
    for tick, cohort in zip(ax.get_xticklabels(), columns):
        tick.set_color(COHORT_COLORS[cohort])
    ax.grid(axis='y', color='#E8E8E8', linewidth=.35)
    ax.set_axisbelow(True)
    return len(selected)


def render_fig1(ctx):
    """Complete genetic input, original RNA geometry, and all47 existing tests."""
    ledger = _frame(ctx, 'magma_all19427',
        ['chromosome', 'official_start', 'official_stop', 'native_P'], 19427)
    groups = _frame(ctx, 'groups47', ['cohort', 'group_family', 'group',
        'n_cells', 'FDR_BH_all47_fixed_MC'], 47)
    atlas220 = _frame(ctx, 'atlas220', ['cell_id', 'UMAP1', 'UMAP2', 'major_lineage'], 75727)
    atlas335 = _frame(ctx, 'atlas335', ['cell_id', 'UMAP1', 'UMAP2', 'major_lineage'], 118426)
    for frame in [atlas220, atlas335]:
        frame['display_lineage'] = frame.major_lineage.astype(str).str.lower()
    fig = plt.figure(figsize=(170/25.4, 193/25.4))
    ax = _axis(fig, [.07, .905, .90, .062]); ax.axis('off')
    labels = [('Published genetic input', 'Fixed MAGMA1000'),
              ('Original single-cell atlases', '47 existing context tests'),
              ('RNA and interpretation', 'Patients, repair, treatments')]
    for i, (top, bottom) in enumerate(labels):
        x = .015 + .337*i
        ax.text(x, .82, top, transform=ax.transAxes, fontsize=6.5, fontweight='bold')
        ax.text(x, .28, bottom, transform=ax.transAxes, fontsize=6.5)
        if i < 2:
            ax.annotate('', xy=(x+.322, .50), xytext=(x+.287, .50),
                        xycoords='axes fraction', arrowprops=dict(arrowstyle='->', lw=.7))
    ctx.panel(ax, 'A', '')
    ax = _axis(fig, [.09, .725, .865, .145])
    p = _numeric(ledger.native_P)
    finite = np.isfinite(p)
    if int(finite.sum()) != 18334 or int((~finite).sum()) != 1093:
        raise ValueError('Preserve full18334 finite /1093 NA ledger')
    if not ((p[finite] > 0) & (p[finite] <= 1)).all():
        raise ValueError('Invalid original gene P value')
    numeric_chr = pd.to_numeric(ledger.chromosome.replace({'X': 23, 'Y': 24}), errors='raise')
    ledger['chr_num'] = numeric_chr.astype(int)
    ledger['midpoint'] = (_numeric(ledger.official_start) + _numeric(ledger.official_stop))/2
    offset = 0.0; ticks = []; ticklabels = []
    for index, chromosome in enumerate(sorted(ledger.chr_num.unique())):
        complete_chr = ledger.loc[ledger.chr_num.eq(chromosome)]
        width = float(complete_chr.official_stop.max())
        points = ledger.loc[ledger.chr_num.eq(chromosome) & finite]
        if len(points):
            ax.scatter(points.midpoint.to_numpy(float)+offset,
                       -np.log10(points.native_P.to_numpy(float)),
                       s=1.2, c=['#7A8FA5', '#313F50'][index % 2],
                       rasterized=True, linewidths=0)
        ticks.append(offset+width/2); ticklabels.append('X' if chromosome == 23 else ('Y' if chromosome == 24 else str(chromosome)))
        offset += width + 2_000_000
    ax.set_xticks(ticks, ticklabels, fontsize=5.5)
    for label in ax.get_xticklabels():
        label.set_rotation(35); label.set_ha('right')
    ax.set_xlabel('Chromosome (gene midpoint, physical spacing within chromosome)', fontsize=6.5, labelpad=1)
    ax.set_ylabel(r'$-\log_{10}$(gene P)', fontsize=6.5, labelpad=1)
    ax.set_xlim(0, offset); ax.set_ylim(bottom=0)
    ctx.panel(ax, 'B', 'Complete published-GWAS gene analysis')
    ax.text(.995, .96, '19,427 genes: 18,334 finite; 1,093 NA',
            transform=ax.transAxes, ha='right', va='top', fontsize=5.5)
    colors = dict(zip(LINEAGE_ORDER, LINEAGE_PALETTE))
    for key, frame, rect, letter, label in [
        ('atlas220', atlas220, [.085, .435, .395, .226], 'C', '220 atlas | 75,727 cells'),
        ('atlas335', atlas335, [.565, .435, .395, .226], 'D', '335 atlas | 118,426 cells')]:
        ax = _axis(fig, rect); _umap(ax, frame, 'display_lineage', colors)
        ctx.panel(ax, letter, label)
    handles = [Line2D([], [], marker='o', linestyle='', color=colors[x],
                      markersize=3, label=LINEAGE_LABELS[x]) for x in LINEAGE_ORDER
               if x in set(atlas220.display_lineage) | set(atlas335.display_lineage)]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .335),
               ncol=5, frameon=False, fontsize=5.5, handletextpad=.35,
               columnspacing=.9, borderaxespad=0, labelspacing=.35)
    cmap = plt.get_cmap('Blues'); norm = Normalize(0, 2.5)
    ax = _axis(fig, [.205, .065, .205, .248])
    present_lineages = set(groups.loc[groups.group_family.eq('lineage'), 'group'].str.lower())
    order = [x for x in LINEAGE_ORDER if x in present_lineages]
    if len(present_lineages.difference(order)):
        raise ValueError('Unregistered lineage, retain source category')
    count = _dot_matrix(ax, groups, 'lineage', order, COHORTS, norm, cmap)
    if count != 32:
        raise ValueError('All32 original lineage tests required')
    ctx.panel(ax, 'E', 'All 32 lineage-context tests')
    cax = fig.add_axes([.465, .076, .011, .220])
    bar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    bar.ax.tick_params(labelsize=5.5, width=.5, length=2)
    bar.set_label(r'$-\log_{10}$(original BH47 q)', fontsize=6.5, labelpad=3)
    ax = _axis(fig, [.73, .065, .18, .248])
    state_order = sorted(groups.loc[groups.group_family.eq('FB_state'), 'group'].unique())
    count = _dot_matrix(ax, groups, 'FB_state', state_order,
                        ['GSE220300', 'GSE335482'], norm, cmap)
    if count != 15:
        raise ValueError('All15 original state tests required')
    if 'FB01' in state_order:
        iy = state_order.index('FB01')
        ax.add_patch(Rectangle((-.48, iy-.48), 1.96, .96, fill=False,
                               edgecolor='#333333', linewidth=.6))
    ctx.panel(ax, 'F', 'All 15 state-context tests')
    fig.text(.07, .017, 'Dot area: log-scaled cells; black outline: q<0.05; ×: original group absent. State IDs are cohort-local.', fontsize=5.5)
    caption = ('Figure 1. Genetic input and complete cellular-context associations. '
        '(A) Public genetic and RNA inputs remain separate. (B) All 18,334 finite native MAGMA gene P values from the full 19,427-gene ledger, '
        'using official gene-midpoint base-pair positions within concatenated chromosomes; 1,093 unavailable entries are retained in source data. '
        'This is gene-level reanalysis of the published GWAS, not discovery of new GWAS loci. (C,D) Every original atlas cell in the previously frozen UMAP, '
        'colored by its preserved lineage; no embedding was refitted or cell score invented. (E,F) All 32 lineage and fifteen state tests retain their original single-family BH47 q values. '
        'Marker area is proportional to 1+log10(cell count); black edges indicate q<0.05. Gray crosses denote groups absent from the original taxonomy, not negative tests. '
        'Display names differ only in capitalization, and Lymphoid and T/NK remain distinct. Cohort-local state IDs do not imply correspondence. '
        'The normal D0 and lesional RNA contexts share one genetic input and do not constitute independent genetic replications.')
    ctx.save(fig, 'Fig1_Genetic_context_and_atlases', caption,
             ['magma_all19427', 'groups47', 'atlas220', 'atlas335'])
    plt.close(fig)


def _state_map(ax, frame, cohort):
    order = sorted(frame.FB_state.astype(str).unique())
    colors = {state: STATE_PALETTE[i % len(STATE_PALETTE)] for i, state in enumerate(order)}
    colors['FB01'] = COHORT_COLORS[cohort]
    _umap(ax, frame, 'FB_state', colors, point_size=.65)
    handles = [Line2D([], [], marker='o', linestyle='', color=colors[x],
               markersize=3, label=x) for x in order]
    ax.legend(handles=handles, frameon=False, ncol=5, fontsize=5.5,
              loc='upper center', bbox_to_anchor=(.5, -.19), handletextpad=.2,
              columnspacing=.55, borderaxespad=0, labelspacing=.25)


def _matrix(ax, values, rows, columns, vmin, vmax, cmap, annotate=False):
    mask = np.ma.masked_invalid(np.asarray(values, dtype=float))
    image = ax.imshow(mask, aspect='auto', cmap=plt.get_cmap(cmap).with_extremes(bad='#BDBDBD'),
                      vmin=vmin, vmax=vmax, interpolation='nearest')
    ax.set_xticks(range(len(columns)), columns, fontsize=5.5)
    ax.set_yticks(range(len(rows)), rows, fontsize=5.5)
    ax.tick_params(length=0)
    if annotate:
        for iy in range(len(rows)):
            for ix in range(len(columns)):
                val = mask[iy, ix]
                if not np.ma.is_masked(val):
                    ax.text(ix, iy, f'{float(val):.2f}', ha='center', va='center',
                            fontsize=5.5, color='white' if val > (vmin+.68*(vmax-vmin)) else '#222222')
    return image


def _patient_matrix(frame, cohort, state_order, patient_order):
    selected = frame.loc[frame.cohort.eq(cohort) & frame.group_family.eq('FB_state') &
                         frame.independent_patient_id.isin(patient_order)].copy()
    if selected[['group', 'independent_patient_id']].duplicated().any():
        raise ValueError('Known patient/state source rows must be unique')
    delta = selected.pivot(index='group', columns='independent_patient_id',
        values='delta_vs_other_FB_same_patient_or_context').reindex(index=state_order, columns=patient_order)
    cells = selected.pivot(index='group', columns='independent_patient_id',
        values='n_cells').reindex(index=state_order, columns=patient_order)
    return np.asarray(delta, dtype=float), np.asarray(cells, dtype=float)


def render_fig2(ctx):
    """Original state identity, all54 correspondence pairs and patient support."""
    fb220 = _frame(ctx, 'fb220', ['cell_id', 'UMAP1', 'UMAP2', 'FB_state'], 26498)
    fb335 = _frame(ctx, 'fb335', ['cell_id', 'UMAP1', 'UMAP2', 'FB_state'], 6774)
    corresponding = _frame(ctx, 'correspondence54', ['state1', 'state2',
        'top50_marker_Jaccard', 'common_gene_centroid_spearman'], 54)
    states = _frame(ctx, 'states15', ['cohort', 'FB_state', 'patient_ids'], 15)
    patients = _frame(ctx, 'patients239', ['cohort', 'group_family', 'group',
        'independent_patient_id', 'n_cells', 'delta_vs_other_FB_same_patient_or_context'], 239)
    external = _frame(ctx, 'external_patient6', ['patient', 'condition', 'n_FB01',
        'n_other_FB_including_Unknown', 'mean_difference_FB01_minus_otherFB'], 6)
    support = _frame(ctx, 'external_support6', ['patient', 'FB01_consensus',
        'other_FB_including_Unknown'], 6)
    if set(external.patient) != set(support.patient):
        raise ValueError('External descriptions and support must identify same six samples')
    support = support.set_index('patient').reindex(external.patient)
    if not np.array_equal(support.FB01_consensus.to_numpy(int), external.n_FB01.to_numpy(int)):
        raise ValueError('Original external state counts disagree')
    if not np.array_equal(support.other_FB_including_Unknown.to_numpy(int), external.n_other_FB_including_Unknown.to_numpy(int)):
        raise ValueError('Original external comparator counts disagree')
    fig = plt.figure(figsize=(170/25.4, 195/25.4))
    for frame, cohort, rect, letter in [
        (fb220, 'GSE220300', [.09, .782, .38, .170], 'A'),
        (fb335, 'GSE335482', [.57, .782, .38, .170], 'B')]:
        ax = _axis(fig, rect); _state_map(ax, frame, cohort)
        ctx.panel(ax, letter, cohort.replace('GSE', '')+' | '+f'{len(frame):,} fibroblasts')
    rows = sorted(corresponding.state1.unique()); columns = sorted(corresponding.state2.unique())
    if len(rows) != 6 or len(columns) != 9:
        raise ValueError('Preserve full6×9 source correspondence')
    for field, rect, cbarrect, letter, title, color in [
        ('top50_marker_Jaccard', [.10, .538, .34, .152], [.457, .548, .008, .132], 'C', 'Top-50-marker Jaccard', 'Blues'),
        ('common_gene_centroid_spearman', [.60, .538, .34, .152], [.957, .548, .008, .132], 'D', 'Common-gene centroid ρ', 'Purples')]:
        ax = _axis(fig, rect)
        table = corresponding.pivot(index='state1', columns='state2', values=field).reindex(index=rows, columns=columns)
        im = _matrix(ax, table.to_numpy(float), rows, columns, 0, 1, color, annotate=True)
        if 'FB01' in rows and 'FB01' in columns:
            ax.add_patch(Rectangle((columns.index('FB01')-.5, rows.index('FB01')-.5), 1, 1,
                                  fill=False, edgecolor='black', linewidth=.9))
        ax.set_xlabel('220 cohort-local state', fontsize=6.5, labelpad=1)
        ax.set_ylabel('335 state', fontsize=6.5, labelpad=1)
        ctx.panel(ax, letter, title)
        cb = fig.colorbar(im, cax=fig.add_axes(cbarrect)); cb.ax.tick_params(labelsize=5.5, length=2)
        cb.set_ticks([0, .5, 1])
    patient_views = []
    for cohort, state_order in [('GSE220300', sorted(fb220.FB_state.unique())),
                                ('GSE335482', sorted(fb335.FB_state.unique()))]:
        source_ids = states.loc[states.cohort.eq(cohort), 'patient_ids'].unique()
        if len(source_ids) != 1:
            raise ValueError('Frozen four-patient order must match across all states')
        patient_order = str(source_ids[0]).split('|')
        if len(patient_order) != 4:
            raise ValueError('Exactly four documented patient slots required')
        delta, celln = _patient_matrix(patients, cohort, state_order, patient_order)
        patient_views.append((cohort, state_order, patient_order, delta, celln))
    vmax = max(float(np.nanmax(np.abs(item[3]))) for item in patient_views)
    if not np.isfinite(vmax) or vmax <= 0:
        raise ValueError('Full source patient contrast range required')
    for item, rect, sidebar, letter in zip(patient_views,
        [[.10, .280, .31, .182], [.60, .280, .31, .182]],
        [[.42, .280, .055, .182], [.92, .280, .055, .182]], ['E', 'F']):
        cohort, state_order, patient_order, delta, celln = item
        ax = _axis(fig, rect)
        im = ax.imshow(np.ma.masked_invalid(delta), aspect='auto',
            cmap=plt.get_cmap('RdBu_r').with_extremes(bad='#BDBDBD'),
            norm=TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), interpolation='nearest')
        ax.set_xticks(range(4), [x.split('::')[-1] for x in patient_order], fontsize=5.5)
        ax.set_yticks(range(len(state_order)), state_order, fontsize=5.5); ax.tick_params(length=0)
        for iy in range(len(state_order)):
            for ix in range(4):
                if np.isfinite(delta[iy, ix]):
                    ax.text(ix, iy, f'{delta[iy, ix]:+.2f}', ha='center', va='center', fontsize=5.5,
                            color='white' if abs(delta[iy, ix]) > .65*vmax else '#222222')
                else:
                    ax.text(ix, iy, 'NA', ha='center', va='center', fontsize=5.5, color='#444444')
        if 'FB01' in state_order:
            iy = state_order.index('FB01'); ax.add_patch(Rectangle((-.5, iy-.5), 4, 1,
                fill=False, edgecolor='black', linewidth=.8))
        ctx.panel(ax, letter, cohort.replace('GSE', '')+' | within-patient Δ')
        n_ax = _axis(fig, sidebar); n_ax.set_xlim(0, 1); n_ax.set_ylim(len(state_order)-.5, -.5); n_ax.axis('off')
        n_ax.text(.5, 1.035, 'cells¹', transform=n_ax.transAxes, ha='center', fontsize=5.5)
        for iy in range(len(state_order)):
            n_ax.text(.5, iy, str(int(np.nansum(celln[iy]))), ha='center', va='center', fontsize=5.5)
    cax = fig.add_axes([.487, .293, .009, .151])
    cb = fig.colorbar(im, cax=cax); cb.ax.tick_params(labelsize=5.5, length=2)
    cb.set_label('State − other FB', fontsize=6.5, labelpad=2)
    ax = _axis(fig, [.13, .066, .30, .141])
    delta = _numeric(external.mean_difference_FB01_minus_otherFB)
    for i, row in enumerate(external.itertuples()):
        if np.isfinite(delta[i]):
            ax.scatter(delta[i], i, marker='o' if row.condition == 'keloid' else 'x',
                       s=17, c=GREEN if row.condition == 'keloid' else GRAY, linewidths=.8)
        else:
            ax.text(.98, i, 'NA (0 state cells)', transform=ax.get_yaxis_transform(),
                    ha='right', va='center', fontsize=5.5, color=GRAY)
    ax.axvline(0, color=GRAY, lw=.6); ax.set_ylim(5.5, -.5)
    ax.set_yticks(range(6), [x.split('::')[-1] for x in external.patient], fontsize=5.5)
    ax.set_xlabel('State − other-FB mean score', fontsize=6.5, labelpad=1)
    ctx.panel(ax, 'G', 'External six-sample descriptions')
    ax = _axis(fig, [.63, .066, .32, .141])
    yy = np.arange(6); n_state = external.n_FB01.to_numpy(int); n_other = external.n_other_FB_including_Unknown.to_numpy(int)
    ax.barh(yy, n_state, color=GREEN, label='FB01-like', height=.65)
    ax.barh(yy, n_other, left=n_state, color='#CCCCCC', label='Other FB', height=.65)
    for iy, count in enumerate(n_state):
        ax.text(float(n_state[iy]+n_other[iy])+28, iy, str(count), va='center', fontsize=5.5, color=GREEN)
    ax.set_ylim(5.5, -.5); ax.set_xlim(0, float(np.max(n_state+n_other))*1.14)
    ax.set_yticks(yy, [x.split('::')[-1] for x in external.patient], fontsize=5.5)
    ax.set_xlabel('Cells (green label = state n)', fontsize=6.5, labelpad=1)
    ax.legend(frameon=False, fontsize=5.5, ncol=2, loc='upper center', bbox_to_anchor=(.5, -.28), borderaxespad=0)
    ctx.panel(ax, 'H', 'External state/comparator coverage')
    fig.text(.07, .015, '¹ Sidebar = total state cells in four documented patients. NA stays missing; no patient-level P values.', fontsize=5.5)
    caption = ('Figure 2. State correspondence and complete patient-level expression descriptions. '
        '(A,B) All previously frozen fibroblast UMAP coordinates colored by their original cohort-local states, with FB01 emphasized. '
        'No state, marker panel or embedding was refitted. (C,D) Every one of the 54 cross-study state pairs retains the original top-50-marker Jaccard '
        'and common-gene centroid Spearman correlation; FB01–FB01 is outlined. Both color scales span zero to one. The broad centroid similarities '
        'do not establish FB01 specificity. (E,F) All fifteen states and all four documented patient slots per cohort retain their original mean state-minus-other-FB '
        'score differences. Gray/NA cells remain unavailable, numbers annotate source differences, and side columns give total state cells across those four patients. '
        'The color range includes every observed contrast. (G,H) All six external samples retain descriptive state contrasts and exact cell coverage, '
        'including ordinary-scar state counts zero, one and one. These sparse observations cannot establish a stable ordinary-scar comparator effect. '
        'There are no cell-bootstrap confidence intervals or new patient-level tests. All RNA sources share the same GWAS and support expression-context transport, not independent genetic replication.')
    ctx.save(fig, 'Fig2_State_identity_and_patient_support', caption,
             ['fb220', 'fb335', 'correspondence54', 'states15', 'patients239',
              'external_patient6', 'external_support6'])
    plt.close(fig)
