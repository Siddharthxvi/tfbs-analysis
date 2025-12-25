#!/usr/bin/env python3
"""
map_fimo_to_promoters.py

Input:
 - promoters.fa : fasta where headers look like:
    >CiarC_t001::NC_011163.1:0-2077
    SEQUENCE...
 - fimo.tsv : standard FIMO TSV with header (tab-separated)
 - tomtom.tsv (optional) : standard TOMTOM output to filter motif IDs (Query_ID column)

Outputs:
 - motif_flanks.fa       (FASTA of motif +/- flank)
 - motif_highlight.html  (HTML with motif spans highlighted inside promoters)
 - motif_summary.csv     (summary table of matches)
 - log messages to stdout about ambiguous/failed mappings

Run using :
python3 highlight_motifs.py --promoters promoters.fa --fimo fimo.tsv --tomtom tomtom.tsv --flank 25 --tomtom_q 0.2 
"""
import argparse
import csv
import html
import sys
from collections import defaultdict

def parse_promoters(fa_path):
    """
    Parse promoters.fa into a dict:
     { accession : [ { 'header': full_header, 'name': gene_name, 'acc': accession,
                       'gen_start': int, 'gen_end': int, 'seq': sequence }, ... ] }
    header format expected: ><gene>::<acc>:<start>-<end>
    """
    promoters = defaultdict(list)
    with open(fa_path) as f:
        header = None
        seq_parts = []
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    # flush previous
                    seq = "".join(seq_parts).replace(" ", "").replace("\r","").replace("\n","")
                    entry = dict(raw_header=header, seq=seq)
                    # parse header
                    # try splitting on '::' and then ':start-end'
                    try:
                        left, right = header[1:].split("::",1)
                        gene_name = left
                        acc_part = right
                    except ValueError:
                        # fallback: try to find accession by the first token after >
                        parts = header[1:].split()
                        gene_name = parts[0]
                        acc_part = parts[0]
                    # extract accession and coords if present
                    acc = None
                    gstart = None
                    gend = None
                    if ":" in acc_part and "-" in acc_part:
                        # e.g. NC_011163.1:0-2077
                        try:
                            acc, coords = acc_part.split(":",1)
                            s,e = coords.split("-",1)
                            gstart = int(s)
                            gend = int(e)
                        except Exception:
                            acc = acc_part
                    else:
                        acc = acc_part
                    entry.update(name=gene_name, acc=acc, gen_start=gstart, gen_end=gend)
                    promoters[entry['acc']].append(entry)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line.strip())
        # flush last
        if header is not None:
            seq = "".join(seq_parts).replace(" ", "").replace("\r","").replace("\n","")
            entry = dict(raw_header=header, seq=seq)
            try:
                left, right = header[1:].split("::",1)
                gene_name = left
                acc_part = right
            except ValueError:
                parts = header[1:].split()
                gene_name = parts[0]
                acc_part = parts[0]
            acc = None
            gstart = None
            gend = None
            if ":" in acc_part and "-" in acc_part:
                try:
                    acc, coords = acc_part.split(":",1)
                    s,e = coords.split("-",1)
                    gstart = int(s)
                    gend = int(e)
                except Exception:
                    acc = acc_part
            else:
                acc = acc_part
            entry.update(name=gene_name, acc=acc, gen_start=gstart, gen_end=gend)
            promoters[entry['acc']].append(entry)
    return promoters

def load_tomtom(tomtom_path, q_thresh=0.2):
    """
    Return set of Query_ID that pass tomtom q-value threshold.
    tomtom.tsv assumed to have header with 'Query_ID' and 'q-value' columns.
    """
    keep = set()
    if not tomtom_path:
        return keep
    with open(tomtom_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if 'Query_ID' not in reader.fieldnames or 'q-value' not in reader.fieldnames:
            # try lowercase or other common names
            pass
        for row in reader:
            try:
                qv = float(row.get('q-value', row.get('q_value','999')))
            except:
                qv = 999.0
            qv = float(qv)
            qid = row.get('Query_ID') or row.get('Query ID') or row.get('Query')
            if qid is None:
                continue
            if qv <= q_thresh:
                keep.add(qid)
    return keep

def map_fimo(promoters, fimo_path, tomtom_keep=None):
    """
    Read fimo.tsv and map to promoters.
    Returns list of match dicts.
    """
    matches = []
    ambiguous = 0
    not_mapped = 0
    with open(fimo_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        # fimo expected header names: motif_id,motif_alt_id,sequence_name,start,stop,strand,score,p-value,q-value,matched_sequence
        for row in reader:
            motif_id = row.get('motif_id') or row.get('motif') or row.get('motif_id')
            seqname = row.get('sequence_name') or row.get('sequence')
            if not motif_id or not seqname:
                continue
            if tomtom_keep is not None and len(tomtom_keep)>0 and motif_id not in tomtom_keep:
                # skip motifs not in tomtom set (user wanted filtering)
                continue
            try:
                fstart = int(row.get('start'))
                fstop  = int(row.get('stop'))
            except Exception:
                # skip malformed line
                continue
            strand = row.get('strand', '.')
            matched_seq = row.get('matched_sequence', '')
            # try to map to promoters by accession match
            promos = promoters.get(seqname)
            if not promos:
                # no promoter for that accession -> not mapped
                not_mapped += 1
                continue
            found = []
            for p in promos:
                pstart = p.get('gen_start')
                pend = p.get('gen_end')
                if pstart is None or pend is None:
                    # can't map by coords if promoter has no genomic coords
                    # fallback: assume this promoter was the sequence_name itself if header acc matches exactly AND lengths equal? skip for now.
                    continue
                # Try two conventions:
                # a) FIMO coordinates are 1-based inclusive (common)
                # b) FIMO coordinates are 0-based (less common)
                # We'll accept hit if motif overlaps promoter region under either convention.
                # Test 1: assume fimo is 1-based inclusive => motif_range = [fstart, fstop]
                if (fstart >= pstart+1) and (fstop <= pend+1):  # promoter stored 0-based start; convert promoter to 1-based
                    # compute relative 1-based inside promoter:
                    rel_start = fstart - pstart   # 1-based
                    rel_end = fstop - pstart
                    found.append((p, rel_start, rel_end, 'fimo_1based'))
                # Test 2: assume fimo is 0-based inclusive
                if (fstart >= pstart) and (fstop <= pend):
                    rel_start = (fstart - pstart) + 1
                    rel_end = (fstop - pstart) + 1
                    found.append((p, rel_start, rel_end, 'fimo_0based'))
            # deduplicate identical findings
            if len(found) == 0:
                not_mapped += 1
                continue
            # keep unique promoter hits (if multiple promoters same accession overlap, record all)
            # but avoid duplicate identical mapping entries
            uniq = []
            seen_keys = set()
            for (p, rs, re, mode) in found:
                key = (p['raw_header'], rs, re)
                if key in seen_keys:
                    continue
                seen_keys.add(key)
                uniq.append((p, rs, re, mode))
            if len(uniq) > 1:
                ambiguous += 1
            for (p, rs, re, mode) in uniq:
                matches.append(dict(
                    motif_id=motif_id,
                    seqname=seqname,
                    fimo_start=fstart,
                    fimo_stop=fstop,
                    strand=strand,
                    matched_sequence=matched_seq,
                    promoter_header=p['raw_header'],
                    promoter_name=p.get('name'),
                    promoter_acc=p.get('acc'),
                    promoter_start=p.get('gen_start'),
                    promoter_end=p.get('gen_end'),
                    rel_start=rs,  # 1-based relative to promoter sequence
                    rel_end=re,
                    mode=mode
                ))
    print(f"Total fimo hits processed: {len(matches)+not_mapped}")
    print(f"Mapped hits: {len(matches)}; Not mapped: {not_mapped}; Ambiguous mappings: {ambiguous}")
    return matches

def write_outputs(matches, promoters, flank, out_prefix="motif"):
    # 1) FASTA of flanks
    fa_out = f"{out_prefix}_flanks.fa"
    html_out = f"{out_prefix}_highlight.html"
    summary_out = f"{out_prefix}_summary.csv"

    # ensure matches grouped by promoter for html highlighting
    promos_hits = defaultdict(list)
    for i,m in enumerate(matches, start=1):
        promos_hits[m['promoter_header']].append((i,m))

    # write FASTA
    with open(fa_out,"w") as fa:
        for i,m in enumerate(matches, start=1):
            hdr = f">{m['promoter_name']}|{m['motif_id']}|{m['promoter_acc']}:{m['fimo_start']}-{m['fimo_stop']}|strand={m['strand']}|rel={m['rel_start']}-{m['rel_end']}"
            # find sequence substring from promoter seq using rel positions (1-based)
            p_seq = None
            # find promoter entry by header
            p_list = promoters.get(m['promoter_acc'],[])
            p_entry = None
            for e in p_list:
                if e['raw_header'] == m['promoter_header']:
                    p_entry = e; break
            if p_entry is None:
                seq_str = m.get('matched_sequence','')
            else:
                p_seq = p_entry['seq']
                # compute indices for slicing (python 0-based)
                # rel_start and rel_end are 1-based inside promoter, convert to 0-based slice:
                start0 = max(0, m['rel_start'] - 1 - flank)
                end0 = min(len(p_seq), m['rel_end'] + flank)  # rel_end is inclusive so end index = rel_end + flank
                seq_str = p_seq[start0:end0]
            fa.write(hdr + "\n")
            # wrap fasta lines at 80
            for pos in range(0, len(seq_str), 80):
                fa.write(seq_str[pos:pos+80] + "\n")

    # write HTML
    with open(html_out,"w") as hf:
        hf.write("<html><head><meta charset='utf-8'><title>Motif highlights</title>\n")
        hf.write("<style>\n.motif { background-color: #ffcccc; font-weight: bold; }\n.promoter { font-family: monospace; white-space: pre; margin-bottom: 1.5em; }\n</style>\n</head><body>\n")
        hf.write("<h1>Motif highlights in promoters</h1>\n")
        for header, hits in promos_hits.items():
            idx, _ = hits[0]
            # get promoter sequence
            p_entry = None
            # find promoter object
            for acc, p_list in promoters.items():
                for e in p_list:
                    if e['raw_header'] == header:
                        p_entry = e; break
                if p_entry: break
            if p_entry is None:
                continue
            p_seq = p_entry['seq']
            # build a list of per-base tags
            tags = [''] * len(p_seq)
            # each hit: place <span> markers
            for i,m in hits:
                start = max(1, m['rel_start'])
                end = max(1, m['rel_end'])
                # convert to 0-based
                s0 = start - 1
                e0 = end  # slice end (exclusive)
                if s0 < 0 or s0 >= len(p_seq):
                    continue
                if e0 > len(p_seq):
                    e0 = len(p_seq)
                # mark with tag id
                for pos in range(s0, e0):
                    tags[pos] = 'X'  # simple marker
            # build highlighted html sequence by grouping
            outline = []
            i_pos = 0
            while i_pos < len(p_seq):
                if tags[i_pos] == 'X':
                    # enter motif span
                    j = i_pos
                    while j < len(p_seq) and tags[j] == 'X':
                        j += 1
                    fragment = html.escape(p_seq[i_pos:j])
                    outline.append(f"<span class='motif'>{fragment}</span>")
                    i_pos = j
                else:
                    j = i_pos
                    while j < len(p_seq) and tags[j] != 'X':
                        j += 1
                    fragment = html.escape(p_seq[i_pos:j])
                    outline.append(fragment)
                    i_pos = j
            hf.write(f"<div class='promoter'><b>{html.escape(header)}</b><br/>{''.join(outline)}</div>\n")
        hf.write("</body></html>\n")

    # write summary CSV
    with open(summary_out,"w",newline='') as sf:
        w = csv.writer(sf)
        w.writerow(['index','motif_id','promoter_header','promoter_name','promoter_acc','promoter_start','promoter_end','fimo_start','fimo_stop','strand','rel_start','rel_end','mode'])
        for i,m in enumerate(matches, start=1):
            w.writerow([i, m['motif_id'], m['promoter_header'], m['promoter_name'], m['promoter_acc'], m['promoter_start'], m['promoter_end'], m['fimo_start'], m['fimo_stop'], m['strand'], m['rel_start'], m['rel_end'], m['mode']])
    print(f"Wrote: {fa_out}, {html_out}, {summary_out}")

def main():
    ap = argparse.ArgumentParser(description="Map FIMO hits to promoter fasta and extract flanks + HTML highlights")
    ap.add_argument("--promoters", required=True, help="promoters.fa")
    ap.add_argument("--fimo", required=True, help="fimo.tsv")
    ap.add_argument("--tomtom", default=None, help="tomtom.tsv (optional) to filter motif IDs)")
    ap.add_argument("--flank", type=int, default=25, help="flank bp upstream and downstream to extract (default 25)")
    ap.add_argument("--tomtom_q", type=float, default=0.2, help="tomtom q-value threshold when filtering (default 0.2)")
    ap.add_argument("--outprefix", default="motif", help="output prefix")
    args = ap.parse_args()

    print("Parsing promoters...")
    promoters = parse_promoters(args.promoters)
    print(f"Parsed promoters for {len(promoters)} accessions. Example accession keys: {', '.join(list(promoters.keys())[:5])}")
    tomtom_keep = set()
    if args.tomtom:
        print(f"Loading tomtom and keeping Query_IDs with q <= {args.tomtom_q} ...")
        tomtom_keep = load_tomtom(args.tomtom, q_thresh=args.tomtom_q)
        print(f"Tomtom kept {len(tomtom_keep)} motif IDs.")
    else:
        tomtom_keep = set()

    print("Mapping fimo hits to promoters...")
    matches = map_fimo(promoters, args.fimo, tomtom_keep=tomtom_keep)
    if len(matches) == 0:
        print("No matches found. Please check that sequence_name in fimo.tsv matches accession in promoters header (e.g. NC_011163.1) and coordinate systems. The script attempts both 0- and 1-based mapping.")
        sys.exit(1)
    write_outputs(matches, promoters, flank=args.flank, out_prefix=args.outprefix)

if __name__ == "__main__":
    main()
