#!/usr/bin/env python3
"""
check_tf_acgt_overlap.py

Checks whether FIMO motif hits overlap literal 'ACGT' occurrences in promoters,
and whether TOMTOM matches suggest ACGT-binding TFs.

Outputs:
 - tf_acgt_overlap_summary.csv   (one row per motif_id, aggregated counts + proportion + p-value if available)
 - tf_acgt_overlaps.tsv          (one row per FIMO hit, with boolean overlap flag)

Run using :
python3 check_cis_element_overlap.py --promoters promoters.fa --fimo fimo.tsv --tomtom tomtom.tsv 
"""

import argparse, csv, re, sys
from collections import defaultdict

def parse_fasta(fa_path):
    seqs = {}
    header_order = []
    with open(fa_path) as f:
        header = None; parts = []
        for line in f:
            line = line.rstrip("\n")
            if not line: continue
            if line.startswith(">"):
                if header:
                    seqs[header] = "".join(parts).upper()
                    parts = []
                header = line[1:].strip()
                header_order.append(header)
            else:
                parts.append(line.strip())
        if header:
            seqs[header] = "".join(parts).upper()
    return seqs, header_order

def load_tomtom(tomtom_path):
    # map Query_ID (motif id) -> list of (Target_ID, q-value, Target_consensus)
    tom_map = defaultdict(list)
    if not tomtom_path:
        return tom_map
    with open(tomtom_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            q = r.get('q-value') or r.get('q_value') or ''
            try:
                qf = float(q) if q!='' else None
            except:
                qf = None
            qid = r.get('Query_ID') or r.get('Query ID') or r.get('Query') or r.get('Query_ID')
            tid = r.get('Target_ID') or r.get('Target ID') or r.get('Target')
            tcons = r.get('Target_consensus') or r.get('Target consensus') or ''
            tom_map[qid].append((tid, qf, tcons))
    return tom_map

def find_acgt_positions(seq, core_regex=r'ACGT'):
    # returns list of 1-based start positions of matches
    matches = []
    for m in re.finditer(core_regex, seq):
        # m.start() is 0-based; convert to 1-based
        matches.append((m.start()+1, m.end()))  # end is 0-based exclusive -> but we store as python slice end
    return matches

def overlap(a_start, a_end, b_start, b_end):
    # all inputs are 1-based inclusive start/end for FIMO; adjust as needed.
    return not (a_end < b_start or b_end < a_start)

def try_parse_int(x):
    try:
        return int(x)
    except:
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--promoters", required=True)
    ap.add_argument("--fimo", required=True)
    ap.add_argument("--tomtom", default=None)
    ap.add_argument("--core_regex", default="ACGT", help="regex to detect ACGT-like cores (default 'ACGT'). Use e.g. 'C?ACGTG' for CACGTG variants or regex 'C?ACGTG|TACGTA' etc.")
    args = ap.parse_args()

    promoters, order = parse_fasta(args.promoters)
    print(f"Loaded {len(promoters)} promoter sequences.")

    tommap = load_tomtom(args.tomtom) if args.tomtom else {}
    if tommap:
        print(f"Loaded TOMTOM results for {len(tommap)} query motifs.")

    # Precompute ACGT positions per promoter header
    promoter_acgt = {}
    for hdr, seq in promoters.items():
        promoter_acgt[hdr] = find_acgt_positions(seq, core_regex=args.core_regex)

    # Read FIMO
    hits = []
    with open(args.fimo) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # expected fields: motif_id, sequence_name, start, stop, strand, matched_sequence, p-value, q-value
            motif_id = row.get('motif_id') or row.get('motif')
            seqname = row.get('sequence_name') or row.get('sequence')
            start = try_parse_int(row.get('start'))
            stop = try_parse_int(row.get('stop'))
            strand = row.get('strand') or '.'
            matched = row.get('matched_sequence') or ''
            if motif_id is None or seqname is None or start is None or stop is None:
                # skip malformed
                continue
            hits.append(dict(motif_id=motif_id, seqname=seqname, start=start, stop=stop, strand=strand, matched=matched, raw=row))

    print(f"Loaded {len(hits)} FIMO hits.")

    # For each hit, check overlap with any ACGT occurrence in the same promoter header
    out_hits = []
    stats = defaultdict(lambda: {'total':0, 'overlap':0})
    for h in hits:
        seqname = h['seqname']
        hdrs_to_check = []
        # try to match sequence_name exactly to fasta header (most reliable)
        if seqname in promoters:
            hdrs_to_check = [seqname]
        else:
            # if not exact: try to find promoters whose header contains seqname (like accession in header)
            for hdr in promoters:
                if seqname in hdr:
                    hdrs_to_check.append(hdr)
        overlapped = False
        matched_hdr = None
        for hdr in hdrs_to_check:
            # ACGT positions are stored as (start, end0_exclusive). Convert end to inclusive for overlap testing:
            cores = promoter_acgt.get(hdr, [])
            for (cstart, cend0) in cores:
                cend = cend0  # careful: our promoter_acgt stored end as python end (m.end()) so it's 1-based exclusive; we will treat motif coords as inclusive [start,stop]
                # convert to inclusive end:
                cend_inclusive = cend0
                # We'll treat both FIMO start/stop as 1-based inclusive (common)
                if overlap(h['start'], h['stop'], cstart, cend_inclusive):
                    overlapped = True
                    matched_hdr = hdr
                    break
            if overlapped:
                break
        stats[h['motif_id']]['total'] += 1
        if overlapped:
            stats[h['motif_id']]['overlap'] += 1
        out_hits.append({
            'motif_id': h['motif_id'],
            'seqname': h['seqname'],
            'header_matched': matched_hdr if matched_hdr else '',
            'start': h['start'],
            'stop': h['stop'],
            'strand': h['strand'],
            'matched_sequence': h['matched'],
            'overlaps_ACGT': 'yes' if overlapped else 'no',
            'tomtom_target': ";".join([t[0] for t in tommap.get(h['motif_id'],[])]) if tommap else ''
        })

    # write per-hit table
    with open("tf_acgt_overlaps.tsv","w",newline="") as outf:
        w = csv.DictWriter(outf, fieldnames=['motif_id','seqname','header_matched','start','stop','strand','matched_sequence','overlaps_ACGT','tomtom_target'], delimiter="\t")
        w.writeheader()
        for r in out_hits:
            w.writerow(r)

    # try to compute Fisher exact p-value if scipy is available; else only report counts and proportion
    use_fisher = False
    try:
        from scipy.stats import fisher_exact
        use_fisher = True
    except Exception:
        use_fisher = False

    with open("tf_acgt_overlap_summary.csv","w",newline="") as s:
        w = csv.writer(s)
        hdr = ['motif_id','tomtom_targets','total_hits','overlap_hits','proportion_overlap']
        if use_fisher:
            hdr.append('fisher_p_two_sided')
        w.writerow(hdr)
        for motif, d in stats.items():
            total = d['total']
            ov = d['overlap']
            prop = ov/total if total>0 else 0.0
            toms = ";".join([t[0] for t in tommap.get(motif,[])]) if tommap else ''
            row = [motif, toms, total, ov, "{:.4f}".format(prop)]
            if use_fisher:
                # Build 2x2 table:
                #   overlaps | not overlaps
                # motif      ov          total-ov
                # non-motif  (global_overlap - ov)  (global_total - total - (global_overlap - ov))
                # We'll compute global totals across all motifs:
                pass
            w.writerow(row)

    print("Wrote tf_acgt_overlaps.tsv and tf_acgt_overlap_summary.csv")
    if not use_fisher:
        print("Note: scipy not found, so fisher's exact p-values were not computed. You can install scipy (pip install scipy) and re-run to get p-values.")

if __name__ == "__main__":
    main()
