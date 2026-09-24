#!/usr/bin/env python
"""Download and construct the paired ClinVar sequence dataset."""
from __future__ import annotations
import argparse, gzip, hashlib, json, re, sys, random
from pathlib import Path
from urllib.parse import urljoin
import requests, pandas as pd
from tqdm import tqdm

CANCER_GENES = sorted(set("TP53 BRCA1 BRCA2 KRAS EGFR PTEN APC PIK3CA BRAF MYC ERBB2 ALK ROS1 MET RET NTRK1 NTRK2 NTRK3 NRAS HRAS CDKN2A RB1 VHL NF1 NF2 STK11 SMAD4 IDH1 IDH2 FLT3 NPM1 KIT JAK2 ABL1 BCR ETV6 RUNX1 PML RARA EZH2 ARID1A ATM CHEK2 PALB2 RAD51C RAD51D MLH1 MSH2 MSH6 PMS2 EPCAM CDH1 CTNNB1 FBXW7 NOTCH1 PTCH1 SUFU WT1".split()))
NON_CANCER_GENES = sorted(set("CFTR HBB PAH GJB2 FBN1 LDLR PKD1 PKD2 MECP2 DMD SMN1 COL1A1 COL1A2 TSC1 TSC2 RYR1 CACNA1S SCN1A SCN5A KCNQ1 KCNH2 OTC GALT HEXA PAH SLC26A4 ATP7B HFE SERPINA1 F8 F9 VWF SDHB SDHD MEN1 RET2 UBE3A FMR1 ASPA LRRK2 GBA1 HTT PMP22 MPZ NF1X".split()))
LABELS = {"Benign":0,"Likely benign":0,"Benign/Likely benign":0,"Benign, Likely benign":0,"Pathogenic":1,"Likely pathogenic":1,"Pathogenic/Likely pathogenic":1,"Pathogenic, Likely pathogenic":1}
BASE = Path(__file__).resolve().parents[2]
GENCODE_URL = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gtf.gz'

def load_gene_strands(gtf: Path) -> dict[str, str]:
    """Return GENCODE gene-symbol -> strand, restricted to gene features."""
    strands = {}
    with gzip.open(gtf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'): continue
            p = line.rstrip().split('\t')
            if len(p) != 9 or p[2] != 'gene': continue
            attrs = dict(re.findall(r'(gene_name|gene_id) "([^"]+)"', p[8]))
            if attrs.get('gene_name') and p[6] in ('+', '-'):
                strands.setdefault(attrs['gene_name'], p[6])
    return strands

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1].upper()

def cap_by_gene(df: pd.DataFrame, cap: int, seed: int = 42) -> pd.DataFrame:
    """Cap each gene while preserving its observed within-gene label ratio."""
    parts=[]
    for gene, group in df.groupby('gene', sort=True):
        if len(group) <= cap:
            parts.append(group)
            continue
        n0=int((group.label==0).sum()); n1=len(group)-n0
        take0=round(cap*n0/len(group)); take1=cap-take0
        g0=group[group.label==0].sample(n=min(take0,n0), random_state=seed)
        g1=group[group.label==1].sample(n=min(take1,n1), random_state=seed)
        parts.append(pd.concat([g0,g1]))
    return pd.concat(parts, ignore_index=True).sample(frac=1, random_state=seed).reset_index(drop=True)

def download(url: str, dest: Path):
    if dest.exists() and dest.stat().st_size > 0: return
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status(); total=int(r.headers.get("content-length",0))
        with open(dest,"wb") as f, tqdm(total=total,unit="B",unit_scale=True,desc=dest.name) as bar:
            for chunk in r.iter_content(1024*1024): f.write(chunk); bar.update(len(chunk))

def discover_clinvar(out: Path, target="2025-07-01") -> tuple[str,str]:
    root="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2025/"
    html=requests.get(root,timeout=60).text
    files=re.findall(r'href="([^"]+clinvar_[0-9]{8}\.vcf\.gz)"',html)
    if not files: files=re.findall(r'(clinvar_[0-9]{8}\.vcf\.gz)',html)
    if not files: raise RuntimeError("No archived ClinVar GRCh38 VCF found")
    files=sorted(set(files)); chosen=min(files,key=lambda x:abs(pd.Timestamp(x[8:16]).date().toordinal()-pd.Timestamp(target).date().toordinal()))
    url=urljoin(root,chosen); download(url,out/chosen); return chosen,url

def parse_info(info: str):
    d={}
    for item in info.split(';'):
        if '=' in item:
            k,v=item.split('=',1); d[k]=v
    return d

def build(vcf: Path, fasta: Path, out: Path, window: int, gtf: Path | None = None, gene_cap: int | None = None):
    import pysam
    genes=CANCER_GENES+NON_CANCER_GENES
    ff=pysam.FastaFile(str(fasta)); rows=[]; stats={"records":0,"excluded_label":0,"gene":0,"ref_fail":0,"unsupported":0}
    strands=load_gene_strands(gtf) if gtf else {}
    # NCBI RefSeq VCF uses chromosome numbers; the downloaded GenBank assembly
    # uses CM accessions for the primary GRCh38 chromosomes.
    chrom_map={str(i):f"CM{663+i-1:06d}.2" for i in range(1,23)}
    chrom_map.update({'X':'CM000685.2','Y':'CM000686.2','MT':'J01415.2','M':'J01415.2'})
    with gzip.open(vcf,"rt") as fh:
        for line in fh:
            if line.startswith('#'): continue
            stats["records"]+=1; p=line.rstrip().split('\t')
            if len(p)<8: continue
            chrom,pos,vid,ref,alts,qual,flt,info=p[:8]; d=parse_info(info); fasta_chrom=chrom_map.get(chrom,chrom)
            geneinfo=d.get("GENEINFO",""); gene_names={x.split(':')[0] for x in geneinfo.split('|') if x}
            gene=next(iter(gene_names & set(genes)),None)
            if not gene: stats["gene"]+=1; continue
            label_raw=d.get("CLNSIG","").replace('_',' ')
            # CLNSIG may contain numeric codes; prefer the human-readable CLNSIG description.
            label_text=d.get("CLNSIG","").replace('%2C',',').replace('|',',').replace('_',' ')
            label=None
            if any(x in label_text.lower() for x in ["pathogenic"]): label=1
            if any(x in label_text.lower() for x in ["benign"]): label=0 if label is None else None
            if label is None: stats["excluded_label"]+=1; continue
            for alt in alts.split(','):
                if len(ref)>50 or len(alt)>50 or any(c not in 'ACGTNacgtn' for c in ref+alt): stats["unsupported"]+=1; continue
                pos0=int(pos)-1; half=window//2; start=max(0,pos0-half); end=start+window
                try: seq=ff.fetch(fasta_chrom,start,end).upper()
                except Exception: stats["ref_fail"]+=1; continue
                rel=pos0-start
                if len(seq)!=window or seq[rel:rel+len(ref)]!=ref.upper(): stats["ref_fail"]+=1; continue
                mut=seq[:rel]+alt.upper()+seq[rel+len(ref):]
                strand=strands.get(gene, '+')
                rows.append(dict(chrom=chrom,pos=int(pos),ref=ref,alt=alt,gene=gene,variant_type="SNV" if len(ref)==len(alt)==1 else "indel",label=label,label_name="Pathogenic" if label else "Benign",genomic_ref_sequence=seq,genomic_mut_sequence=mut,ref_sequence=reverse_complement(seq) if strand=='-' else seq,mut_sequence=reverse_complement(mut) if strand=='-' else mut,strand=strand,source_release=vcf.name))
    columns=['chrom','pos','ref','alt','gene','variant_type','label','label_name','genomic_ref_sequence','genomic_mut_sequence','ref_sequence','mut_sequence','strand','source_release']
    df=pd.DataFrame(rows,columns=columns).drop_duplicates(["chrom","pos","ref","alt"])
    uncapped_rows=len(df)
    if gene_cap:
        df=cap_by_gene(df, gene_cap)
    out.parent.mkdir(parents=True,exist_ok=True); df.to_parquet(out,index=False)
    stats.update(uncapped_rows=uncapped_rows, gene_cap=gene_cap, final_rows=len(df),benign=int((df['label']==0).sum()),pathogenic=int((df['label']==1).sum()),strand_counts=df['strand'].value_counts().to_dict())
    print(json.dumps({"window":window,"stats":stats,"anchor":"validation anchor is 8,161 (4,818 benign / 3,343 pathogenic); do not force-match"},indent=2))
    return df

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--window',type=int,default=100); ap.add_argument('--vcf'); ap.add_argument('--fasta'); ap.add_argument('--gtf'); ap.add_argument('--gene-cap',type=int); ap.add_argument('--out'); args=ap.parse_args()
    raw=BASE/'data/raw'; raw.mkdir(parents=True,exist_ok=True)
    fasta=Path(args.fasta) if args.fasta else raw/'GRCh38_no_alt_genomic.fna.gz'
    if not fasta.exists(): download('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.29_GRCh38.p14/GCA_000001405.29_GRCh38.p14_genomic.fna.gz',fasta)
    # pysam needs an indexed FASTA; decompress to a stable local path and index it.
    fa=raw/'GRCh38_no_alt_genomic.fna'
    if not fa.exists():
        with gzip.open(fasta,'rb') as src, open(fa,'wb') as dst:
            for chunk in iter(lambda:src.read(1024*1024),b''): dst.write(chunk)
    import pysam
    if not Path(str(fa)+'.fai').exists(): pysam.faidx(str(fa))
    gtf=Path(args.gtf) if args.gtf else raw/'gencode.v46.annotation.gtf.gz'
    if not gtf.exists(): download(GENCODE_URL,gtf)
    vcf=Path(args.vcf) if args.vcf else raw/'clinvar_selected.vcf.gz'
    if not vcf.exists():
        name,url=discover_clinvar(raw); vcf=raw/name
        print(json.dumps({"clinvar":name,"url":url}))
    build(vcf,fa,Path(args.out) if args.out else BASE/f'data/processed/variants_{args.window}bp.parquet',args.window,gtf,args.gene_cap)
if __name__=='__main__': main()
