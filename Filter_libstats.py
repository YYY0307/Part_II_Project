#!/usr/bin/env python3
"""
Usage:
Filter TE stats TSV from libstats.R with class-specific criteria.

Universal:
  - maxTR < 33% of TE_len

List of family is compiled by extracting families in the libraries using Onix function grapping text after # Helper
Run #Helper first and modify the family name to future data set if applicable

Class-specific rules:

  LTR retrotransposons (must have LTR structure):
    Families (from name): LTR, CLASSI/LTR, LTR/BELPAO, LTR/COPIA, LTR/Copia,
                          LTR/ERVK, LTR/GYPSY, LTR/Gypsy,
                          LTR/LARD, LTR/Ngaro, LTR/Pao, LTR/TRIM
    Conditions:
      - struct_type == 'LTR'
      - left_gap + right_gap < 15
      - 3–15 kb (TE_len)
      - struct_len = maxTR
      - orf1 > 2500bp

  TIR DNA transposons (must have TIR structure):
    Families (from name; canonical TIR groups, RM2-style and CLASSII-style):
      DNA, DNA/CMC*, DNA/CMC-*, DNA/CMC-EnSpm, DNA/CMC-Transib
      DNA/HAT, DNA/hAT*, DNA/MULE*, DNA/Merlin,
      DNA/P, DNA/PIF-*, DNA/PIFHARBINGER,
      DNA/PIGGYBAC, DNA/PiggyBac,
      DNA/TC1MARINER, DNA/TRANSIB,
      DNA/TcMar-*, DNA/Zisupton, DNA/MITE
      CLASSII/DNA, CLASSII/DNA/CMC*, CLASSII/DNA/HAT, CLASSII/DNA/hAT*,
      CLASSII/DNA/MULE*, CLASSII/DNA/Merlin,
      CLASSII/DNA/P, CLASSII/DNA/PIF-*, CLASSII/DNA/PIFHARBINGER,
      CLASSII/DNA/PIGGYBAC, CLASSII/DNA/PiggyBac,
      CLASSII/DNA/TC1MARINER, CLASSII/DNA/TRANSIB,
      CLASSII/DNA/TcMar-*, CLASSII/DNA/Zisupton, CLASSII/MITE 
    Conditions:
      - struct_type == 'TIR'
      - left_gap + right_gap < 15


  Mavericks:
    Families: DNA/Maverick, MAVERICK
    Conditions:
      - struct_type == 'TIR'
      - left_gap + right_gap < 15
      - 12–25 kb (TE_len)

  LINEs:
    Families: LINE, LINE/*
    Conditions:
      - 3.5–12.5 kb (TE_len)
      - orf1 > 2000

  Helitrons and Helentrons:
    Families: *Helitron* or *Helentron* (substring match)
    Conditions:
      - TE_len > 9 kb
      - orf1 > 2.5 kb

  MITEs, SINEs, CRYPTON, PLE, Unknown, rRNA:
    - SINE, CRYPTON, PLE, Unknown, rRNA only have universal filter.

Usage:
  python Filter_libstats.py input.stats [output_pass.stats]
Creates:
  <output_pass> and a corresponding <output_fail>.
"""

import sys
import pandas as pd
import numpy as np
import argparse


def main():
    parser = argparse.ArgumentParser(description='Advanced TE stats filter from libstats.R')
    parser.add_argument('input', help='Input .stats TSV')
    parser.add_argument('output', nargs='?', default=None,
                        help='Output TSV for entries that PASS filters')
    parser.add_argument('--fasta', help='Original FASTA to extract filtered seqs for PASS set')
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep='\t')

    numeric_cols = ['TE_len', 'maxTR', 'struct_len', 'left_gap', 'right_gap', 'orf1', 'orf2', 'orf3']
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    
    # Helper: extract Family from name (RepeatMasker Binary Classification expected)
    def extract_family(name):
        if pd.isna(name):
            return ''
        s = str(name)
        if '#' in s:
            fam = s.split('#', 1)[1]
        else:
            fam = s
        return fam

    df['Family'] = df['name'].apply(extract_family)

    
    # Universal filter: maxTR < 33% TE_len
    
    base_mask = df['maxTR'] < 0.33 * df['TE_len']


    # Class-specific filters:
    
    class_mask = pd.Series(True, index=df.index)

    # 1) LTR retrotransposons (by Family)
    ltr_family_flag = df['Family'].str.contains(
        r'^LTR\b|^LTR/|^CLASSI/LTR\b|^CLASSI/LTR/|/BELPAO|/COPIA|/Copia|/GYPSY|/Gypsy|/ERVK|/LARD|/Ngaro|/Pao|/TRIM',
        case=False, regex=True, na=False
    )

    ltr_mask = (
        ltr_family_flag &
        (df['struct_type'] == 'LTR') &
        (df['left_gap'] + df['right_gap'] < 15) &
        (df['TE_len'] >= 3000) & (df['TE_len'] <= 15000) &
        np.isclose(df['struct_len'], df['maxTR'], atol=10) &
        (df['orf1'] > 2500)
    )

    # An LTR retrotransposon must have struct_type LTR
    class_mask &= ~ltr_family_flag | ltr_mask

    # 2) TIR DNA transposons (by Family; require struct_type == 'TIR')
    tir_family_flag = df['Family'].str.contains(
        r'^DNA\b|^DNA/CMC\b|^DNA/CMC-|^DNA/HAT\b|^DNA/hAT\b|^DNA/hAT-|^DNA/MULE\b|^DNA/MULE-|'
        r'^DNA/Merlin\b|^DNA/P\b|^DNA/PIF-|^DNA/PIFHARBINGER\b|'
        r'^DNA/PIGGYBAC\b|^DNA/PiggyBac\b|'
        r'^DNA/TC1MARINER\b|^DNA/TRANSIB\b|^DNA/TcMar-|^DNA/Zisupton\b|^DNA/MITE\b|'
        r'^CLASSII/DNA\b|^CLASSII/DNA/CMC\b|^CLASSII/DNA/CMC-|^CLASSII/DNA/HAT\b|'
        r'^CLASSII/DNA/hAT\b|^CLASSII/DNA/hAT-|^CLASSII/DNA/MULE\b|^CLASSII/DNA/MULE-|'
        r'^CLASSII/DNA/Merlin\b|^CLASSII/DNA/P\b|^CLASSII/DNA/PIF-|^CLASSII/DNA/PIFHARBINGER\b|'
        r'^CLASSII/DNA/PIGGYBAC\b|^CLASSII/DNA/PiggyBac\b|'
        r'^CLASSII/DNA/TC1MARINER\b|^CLASSII/DNA/TRANSIB\b|^CLASSII/DNA/TcMar-|'
        r'^CLASSII/DNA/Zisupton\b|^CLASSII/MITE\b',
        case=False, regex=True, na=False
    )

    tir_mask = (
        tir_family_flag &
        (df['struct_type'] == 'TIR') &
        ((df['left_gap'] + df['right_gap']) < 15) &
        np.isclose(df['struct_len'], df['maxTR'], atol=10)

    )

    # If in DNA family → must have structure TIR;
    class_mask &= ~tir_family_flag | tir_mask


    # 3) Mavericks: must be TIR, gaps < 15, struct_len 12–25 kb
    mav_family_flag = df['Family'].str.contains(
        r'Maverick', case=False, regex=True, na=False
    )

    maverick_mask = (
        mav_family_flag &
        (df['struct_type'] == 'TIR') &
        (df['left_gap'] + df['right_gap'] < 15) &
        (df['TE_len'] >= 12000) & (df['TE_len'] <= 25000)
    )

    class_mask &= ~mav_family_flag | maverick_mask

    # 4) LINEs: 3.5–12.5 kb, ORF1 > 2000
    line_family_flag = df['Family'].str.contains(
        r'^LINE\b|^LINE/', case=False, regex=True, na=False
    )

    line_mask = (
        line_family_flag &
        (df['TE_len'] >= 3500) & (df['TE_len'] <= 12500) &
        (df['orf1'] > 2000)
    )

    class_mask &= ~line_family_flag | line_mask

    # 5) Helitrons and Helentrons: TE_len > 9 kb, ORF1 > 2.5 kb
    helitron_family_flag = df['Family'].str.contains(
        r'Helitron|Helentron', case=False, regex=True, na=False
    )

    helitron_mask = (
        helitron_family_flag &
        (df['TE_len'] > 9000) &
        (df['orf1'] > 2500)
    )

    class_mask &= ~helitron_family_flag | helitron_mask

    # 6) Unknown: always fail
    unknown_flag = df['Family'].str.startswith('Unknown', na=False)
    class_mask &= ~unknown_flag


    # Final mask: universal AND class-specific

    full_mask = base_mask & class_mask

    passed = df[full_mask].dropna(subset=['TE_len', 'maxTR'])
    failed = df[~full_mask].copy()

   
    # Output filenames
 
    if args.output is None:
        out_pass = args.input.replace('.stats', '_filtered_pass.stats')
    else:
        out_pass = args.output

    if '_pass' in out_pass:
        out_fail = out_pass.replace('_pass', '_fail')
    else:
        out_fail = out_pass.replace('.stats', '_fail.stats')

    passed.to_csv(out_pass, sep='\t', index=False)
    failed.to_csv(out_fail, sep='\t', index=False)

 
    # Extract FASTA for passed 
   
    if args.fasta:
        names = passed['name'].astype(str).str.replace('^>', '', regex=True).tolist()
        import subprocess, tempfile, os
        with tempfile.NamedTemporaryFile('w', delete=False) as tmp:
            for n in names:
                tmp.write(n + '\n')
            tmp_name = tmp.name
        fasta_out = args.fasta.replace('.fa', '_filtered_pass.fa')
        cmd = f"seqtk subseq {args.fasta} {tmp_name} > {fasta_out}"
        subprocess.run(cmd, shell=True, check=True)
        import os
        os.unlink(tmp_name)
        print(f"Filtered FASTA (pass): {fasta_out}")

    
    # Summary (output in Terminak)
  
    print(f"Total:  {len(df)}")
    print(f"Passed: {len(passed)} → {out_pass}")
    print(f"Failed: {len(failed)} → {out_fail}")
    print("\nPer-family breakdown (pass):")
    print(passed['Family'].value_counts())


if __name__ == '__main__':
    main()
