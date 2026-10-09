#!/usr/bin/env python3
"""Validate an original-compatible run without changing supplied DNA conventions."""
import argparse, csv, hashlib, json, re
from pathlib import Path

def sha(path):
    h=hashlib.sha256()
    with open(str(path),'rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def rc(sequence):
    return sequence.translate(str.maketrans('ACGTN','TGCAN'))[::-1]

def original_linker_common(sequence, min_length=15):
    """Match upstream linker_common.R, called by read_sample_files.R."""
    parts=re.split('N+', sequence)
    while parts and parts[-1]=='':parts.pop() # R strsplit drops trailing empty parts
    if len(parts)==2:
        common=parts[1]
    elif len(parts)==1:
        if len(parts[0])<=min_length:raise ValueError('Linker must be longer than 15 bases')
        common=parts[0][-min_length:]
    else:raise ValueError('Unsupported linker N layout')
    if len(common)<min_length:raise ValueError('Linker suffix must contain at least 15 bases')
    return rc(common)

def read_sheet(path):
    lines=Path(path).read_text().splitlines()
    start=lines.index('[Data]')+1 if '[Data]' in lines else 0
    return list(csv.DictReader(lines[start:],delimiter='\t' if '\t' in lines[start] else ','))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('samplesheet');p.add_argument('runfolder');p.add_argument('references')
    p.add_argument('source_root');p.add_argument('--fastq-manifest',default='')
    p.add_argument('--annotation',choices=['bushman','off'],default='bushman')
    p.add_argument('--demux',choices=['original','fixed','golay'],default='original')
    p.add_argument('--input-type',choices=['BCL','FASTQ'],default='FASTQ')
    p.add_argument('--runtime-image',required=True)
    p.add_argument('--blat-binary',required=True)
    a=p.parse_args()
    rows=read_sheet(a.samplesheet)
    if not rows:raise ValueError('Empty samplesheet')
    refs=json.loads(Path(a.references).read_text())
    shortest=min(len(r.get('index','')) for r in rows)
    samples=set(); indices=set();files={};facts={}
    if a.fastq_manifest:
        for item in read_sheet(a.fastq_manifest):
            sid=item.get('Sample_ID') or item.get('sample')
            if sid in files:raise ValueError('Duplicate FASTQ sample: '+sid)
            files[sid]={mate:str(Path(item[mate]).resolve()) for mate in ['r1','r2']}
            if a.demux=='original':
                if not item.get('index_fastq'):
                    raise ValueError('Original processing requires index_fastq in the assigned FASTQ manifest; R1/R2 alone cannot reproduce index-quality filtering and Golay correction')
                files[sid]['index_fastq']=str(Path(item['index_fastq']).resolve())
            if not all(Path(x).is_file() for x in files[sid].values()):raise ValueError('Missing FASTQ for '+sid)
    Path('samples').mkdir()
    normalized=[]
    for original in rows:
        row=dict(original)
        row['Sample_ID']=row.get('Sample_ID') or row.get('alias','')
        sid=row['Sample_ID']
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',sid) or sid in samples:raise ValueError('Invalid/duplicate sample')
        samples.add(sid)
        row['ltrbit']=row.get('ltrbit') or row.get('ltrBit','')
        row['Sample_Project']=row.get('Sample_Project') or 'INSPIIRED'
        row['Specimen_ID']=row.get('Specimen_ID') or row.get('GTSP') or sid
        row['Replicate_ID']=row.get('Replicate_ID') or sid
        row['Report_Group']=row.get('Report_Group') or (
            row.get('Trial','trial')+'.'+row['Patient'] if row.get('Patient') else row['Sample_Project'])
        for key in ['Specimen_ID','Replicate_ID','Report_Group']:
            if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',row[key]):raise ValueError('Invalid '+key)
        for key in ['mingDNA','maxFragLength']:
            if not row.get(key,'').isdigit() or int(row[key])<1:raise ValueError('Invalid '+key)
        if not row.get('maxAlignStart','').isdigit():raise ValueError('Invalid maxAlignStart')
        if not 0<=float(row['minPctIdent'])<=100:raise ValueError('Invalid identity threshold')
        if not row.get('linkerSequence'):
            if len(row['index'])-shortest not in [0,1]:raise ValueError('Unsupported inline barcode length difference')
            row['linkerSequence']='N'*(12+len(row['index'])-shortest)+row['common_linker']
        row['requested_linkerCommon']=row.get('linkerCommon','')
        if a.demux in ('original','golay'):
            # Upstream metadata processing overwrites even explicitly supplied
            # linkerCommon. Compare effective metadata, not the raw CSV column.
            row['linkerCommon']=original_linker_common(row['linkerSequence'])
            row['linkerCommon_source']='upstream linker_common(linkerSequence)'
        else:
            row['linkerCommon']=row.get('linkerCommon') or rc(row['common_linker'])
            row['linkerCommon_source']='fixed adapter explicit value or reverse complement'
        for key,default in [('qualityThreshold','10'),('badQualityBases','5'),('qualitySlidingWindow','10')]:
            row[key]=row.get(key) or default
            if float(row[key])<0:raise ValueError('Invalid '+key)
        if not row['qualitySlidingWindow'].isdigit() or int(row['qualitySlidingWindow'])<1:
            raise ValueError('Invalid qualitySlidingWindow')
        for key in ['primer','ltrbit','largeLTRFrag','linkerSequence','linkerCommon']:
            if not re.fullmatch('[ACGTN]+',row[key]):raise ValueError('Invalid DNA field '+key)
        if 'N' not in row['linkerSequence']:raise ValueError('linkerSequence must contain the original primerID N region')
        vector=Path(row['vectorSeq'])
        if not vector.is_absolute():vector=Path(a.runfolder)/vector
        row['vector_fasta']=str(vector.resolve())
        if not vector.is_file():raise ValueError('Missing vector '+sid)
        genome=row['refGenome']
        if genome not in refs or not refs[genome].get('twobit'):
            raise ValueError('Reference manifest must explicitly select the original-compatible twoBit for '+genome)
        entry=refs[genome]
        path=Path(entry['twobit'])
        if not path.is_absolute() or not path.is_file():raise ValueError('twoBit must be an existing absolute path')
        row['reference_twobit']=str(path.resolve())
        row['reference_seqinfo_rds']=entry.get('seqinfo_rds','')
        if row['reference_seqinfo_rds']:
            seqinfo=Path(row['reference_seqinfo_rds'])
            if not seqinfo.is_absolute() or not seqinfo.is_file():raise ValueError('seqinfo_rds must be an existing absolute path')
            actual=sha(seqinfo)
            if entry.get('seqinfo_sha256') and actual!=entry['seqinfo_sha256']:raise ValueError('Seqinfo SHA256 mismatch')
            facts[str(seqinfo)]=dict(sha256=actual,kind='frozen original reference genome/circularity labels')
        if str(path) not in facts:
            actual=sha(path)
            if entry.get('sha256') and actual!=entry['sha256']:raise ValueError('Reference SHA256 mismatch')
            facts[str(path)]=dict(sha256=actual,assembly=genome,sequence_policy=entry.get('sequence_policy','unspecified'))
        row['refseq_rds']=entry.get('refseq_rds','')
        if a.annotation=='bushman' and not Path(row['refseq_rds']).is_file():
            raise ValueError('Provide the frozen assembly-matched RefSeq GRanges RDS, or explicitly select --annotation_mode off')
        if row['refseq_rds']:facts[row['refseq_rds']]=dict(sha256=sha(row['refseq_rds']),kind='RefSeq annotation')
        if files:row.update(files[sid])
        if (a.demux=='fixed' or (a.demux=='original' and a.input_type=='BCL')) and not files:
            key=(row['index'][:shortest],row['index2'])
            if key in indices:raise ValueError('Normalized index collision')
            indices.add(key)
        if a.demux in ('golay','original'):
            barcode=row.get('bcSeq') or row.get('index2','')
            if not re.fullmatch('[ACGT]{12}',barcode) or barcode in indices:
                raise ValueError('Golay barcodes must be unique 12-base sequences')
            row['bcSeq']=barcode
            indices.add(barcode)
        Path('samples/'+sid+'.json').write_text(json.dumps(row,indent=2)+'\n')
        normalized.append(row)
    if files and set(files)!=samples:raise ValueError('FASTQ manifest must contain exactly the expected samples')
    fields=list(dict.fromkeys(key for row in normalized for key in row))
    with open('validated_samples.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t');w.writeheader();w.writerows(normalized)
    root=Path(a.source_root)
    frozen=json.loads((root/'vendor/bushman/PROVENANCE.json').read_text())
    if any(sha(root/'vendor/bushman'/item['component']/item['file'])!=item['sha256'] for item in frozen['files']):
        raise ValueError('Frozen upstream scientific source has changed')
    source_files=[f for folder in ['vendor/bushman']
                  for f in (root/folder).rglob('*') if f.is_file()]
    source_files += [root/'bin/bushman_inputs.py',root/'vendor/bushman_r_sources/runtime_lock.json']
    identity=json.loads((root/'pipeline_identity.json').read_text())
    provenance=dict(analysis=identity['pipeline']+': original-compatible INSPIIRED',pipeline_version=identity['version'],sources={str(f.relative_to(root)):sha(f) for f in source_files},
        references=facts,vectors={r['vector_fasta']:sha(r['vector_fasta']) for r in normalized},
        samplesheet_sha256=sha(a.samplesheet),input_adapter=('assigned FASTQ followed by original Golay' if a.demux=='original' else 'per-sample FASTQ') if files else a.demux,
        original_index_processing=(a.demux in ('original','golay')),
        annotation=a.annotation,marker_convention='original/golay linkerCommon derived from linkerSequence as in upstream effective metadata; supplied value audited as requested_linkerCommon')
    provenance['report_scope']='Report_Group (Trial.Patient when supplied; otherwise Sample_Project) and reference, before specimen splitting'
    provenance['source_scope']='Frozen Bushman and input-validation dependencies; active implementation fingerprints are in pipeline_implementation.json'
    environment=[Path(a.runtime_image),Path(a.blat_binary),root/'vendor/bushman_r_sources/runtime_lock.json']
    if not all(path.is_file() for path in environment):raise ValueError('Missing runtime image, BLAT binary or dependency lock')
    provenance['runtime']={str(path.resolve()):sha(path) for path in environment}
    Path('input_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')

if __name__=='__main__':main()
