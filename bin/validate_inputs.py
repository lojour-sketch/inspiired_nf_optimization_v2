#!/usr/bin/env python3
"""Validate library metadata and fingerprint exact reference/vector inputs."""
import argparse
import csv
import glob
import hashlib
import json
import re
import sys
from pathlib import Path

TXDB = {'hg38':'TxDb.Hsapiens.UCSC.hg38.refGene', 'hg19':'TxDb.Hsapiens.UCSC.hg19.knownGene',
        'hg18':'TxDb.Hsapiens.UCSC.hg18.knownGene'}
REQUIRED = ['Sample_ID','index','index2','common_linker','primer','ltrbit','largeLTRFrag',
            'Sample_Project','mingDNA','minPctIdent','maxAlignStart','maxFragLength','refGenome','vectorSeq']


def fingerprint(path):
    digest = hashlib.sha256()
    with open(path,'rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''): digest.update(block)
    return digest.hexdigest()


def read_samplesheet(path):
    lines = Path(path).read_text().splitlines()
    start = lines.index('[Data]') + 1 if '[Data]' in lines else 0
    rows = list(csv.DictReader(lines[start:]))
    if not rows or not set(REQUIRED).issubset(rows[0]): raise ValueError('Missing samplesheet rows/required columns')
    return rows


def validate_rows(rows):
    ids, indices = set(), set()
    if max(len(row['index']) for row in rows)-min(len(row['index']) for row in rows)>1:
        raise ValueError('Only one-base index normalization (12/13 nt UMI layout) is supported')
    shortest_index = min(len(row['index']) for row in rows)
    for row in rows:
        sample = row['Sample_ID']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',sample) or sample in ids:
            raise ValueError('Invalid/duplicate Sample_ID: '+sample)
        ids.add(sample)
        index = (row['index'][:shortest_index], row['index2'])
        if index in indices: raise ValueError('Duplicate index pair after index-length normalization: '+sample)
        indices.add(index)
        for column in ['index','index2','common_linker','primer','ltrbit','largeLTRFrag']:
            if not re.fullmatch('[ACGTN]+',row[column]): raise ValueError('Invalid DNA field '+column+' for '+sample)
        for column in ['mingDNA','maxAlignStart','maxFragLength']:
            if not row[column].isdigit() or int(row[column]) <= 0: raise ValueError('Invalid '+column+' for '+sample)
        if not 0 <= float(row['minPctIdent']) < 100: raise ValueError('Identity must be >=0 and <100: '+sample)
        if row['refGenome'] not in TXDB: raise ValueError('Unsupported reference assembly: '+row['refGenome'])
        for key in ['Specimen_ID','Replicate_ID']:
            row[key] = row.get(key) or sample
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*',row[key]): raise ValueError('Invalid '+key)
    if len({row['Sample_Project'] for row in rows}) != 1:
        raise ValueError('One Sample_Project per run is required by BCL project routing')
    return rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('samplesheet'); parser.add_argument('runfolder')
    parser.add_argument('--reference-manifest', default='')
    parser.add_argument('--project-root', default='')
    parser.add_argument('--runtime-config', default='')
    args=parser.parse_args()
    rows=validate_rows(read_samplesheet(args.samplesheet))
    references=json.loads(Path(args.reference_manifest).read_text()) if args.reference_manifest else {}
    facts, warnings = {}, []
    for row in rows:
        genome=row['refGenome']
        entry=references.get(genome,{})
        paths = [entry['fasta']] if entry.get('fasta') else glob.glob(str(Path(args.runfolder).parent/(genome+'*.fa')))
        if entry.get('fasta') and not Path(paths[0]).is_absolute():
            raise ValueError('Reference manifest FASTA paths must be absolute: '+genome)
        if len(paths)!=1 or not Path(paths[0]).is_file(): raise ValueError('Select exactly one existing FASTA for '+genome)
        fasta=Path(paths[0]).resolve()
        row['reference_fasta']=str(fasta)
        row['reference_txdb']=entry.get('txdb',TXDB[genome])
        if row['reference_txdb'] != TXDB[genome]: raise ValueError('TxDb/assembly mismatch for '+genome)
        if str(fasta) not in facts:
            dictionary=[]
            with fasta.open() as handle:
                name, length=None,0
                for line in handle:
                    if line.startswith('>'):
                        if name is not None: dictionary.append({'name':name,'length':length})
                        name=line[1:].split()[0]; length=0
                    else: length+=len(line.strip())
                if name is not None: dictionary.append({'name':name,'length':length})
            if not dictionary or any(x['length'] == 0 for x in dictionary) or len({x['name'] for x in dictionary}) != len(dictionary):
                raise ValueError('Empty FASTA/duplicate sequence names')
            actual=fingerprint(fasta)
            if entry.get('sha256') and entry['sha256'] != actual: raise ValueError('Reference checksum mismatch')
            facts[str(fasta)]=dict(sha256=actual,assembly=genome,accession=entry.get('accession','unspecified'),
                sequence_policy=entry.get('sequence_policy','unspecified'),dictionary=dictionary)
            if not entry: warnings.append(genome+': inferred legacy FASTA; sequence set is not harmonized automatically')
            if genome=='hg38' and len(dictionary)<=25: warnings.append('hg38 primary-only reference: alternate/patch ambiguity is unassessed')
        vector=Path(row['vectorSeq'])
        if not vector.is_absolute(): vector=Path(args.runfolder)/vector
        if not vector.is_file(): raise ValueError('Missing vector FASTA for '+row['Sample_ID'])
        row['vector_fasta']=str(vector.resolve())
    fields=REQUIRED+['Specimen_ID','Replicate_ID','reference_fasta','reference_txdb','vector_fasta']
    with open('validated_samples.tsv','w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
    provenance = dict(schema_version=3, python_version=sys.version,
        samplesheet_sha256=fingerprint(args.samplesheet),references=facts,
        vectors={row['vector_fasta']:fingerprint(row['vector_fasta']) for row in rows},warnings=warnings)
    if args.project_root:
        root = Path(args.project_root)
        provenance['source_sha256'] = {str(p.relative_to(root)): fingerprint(p)
            for directory in ['bin','modules','subworkflows','conf']
            for p in sorted((root/directory).rglob('*')) if p.is_file() and p.suffix in ['.nf','.py','.R','.config','.json']}
        for name in ['main.nf','nextflow.config']:
            provenance['source_sha256'][name] = fingerprint(root/name)
        provenance['local_container_sha256'] = {str(p.relative_to(root)): fingerprint(p)
            for p in sorted((root/'containers').glob('*.sif'))}
    if args.runtime_config:
        provenance['runtime'] = json.loads(Path(args.runtime_config).read_text())
    Path('input_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')


if __name__=='__main__': main()
