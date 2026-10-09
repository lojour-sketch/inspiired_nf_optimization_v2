#!/usr/bin/env python3
"""Add real Golay indexes and native/multi-project entry fixtures to BLAT truth reads."""
import csv, gzip, sys, types
from pathlib import Path

fixture=Path(sys.argv[1]).resolve();pipeline=Path(sys.argv[2]).resolve()
source=pipeline/'vendor/bushman/intSiteCaller/errorCorrectIndices/golay.py'
code=source.read_text().replace('numpy.array(map(int,bitstring))','numpy.array(list(map(int,bitstring)))')
golay=types.ModuleType('golay');exec(compile(code,str(source),'exec'),golay.__dict__)
rows=list(csv.DictReader((fixture/'samples.csv').open()))
manifest=list(csv.DictReader((fixture/'fastqs.csv').open()))
raw=fixture/'raw';raw.mkdir(exist_ok=True)
streams={mate:gzip.open(str(raw/('Undetermined_S0_'+mate+'_001.fastq.gz')),'wt') for mate in ['R1','R2','I2']}
for number,row in enumerate(rows,1):
    sid=row['Sample_ID'];barcode=golay.encode([int(x) for x in format(number,'012b')])
    row.update(index2=barcode,bcSeq=barcode,linkerSequence='N'*12+row['common_linker'],
               linkerCommon=row['common_linker'],qualityThreshold='10',badQualityBases='5',qualitySlidingWindow='10',Report_Group='truth')
    index=fixture/(sid+'.I2.fq.gz');wrong_index=fixture/(sid+'.I1.fq.gz')
    with gzip.open(str(manifest[number-1]['r1']),'rt') as first, gzip.open(str(manifest[number-1]['r2']),'rt') as second, \
         gzip.open(str(index),'wt') as target,gzip.open(str(wrong_index),'wt') as wrong:
        while True:
            r1=[first.readline() for _ in range(4)]
            if not r1[0]:break
            r2=[second.readline() for _ in range(4)];assert r1[0]==r2[0]
            ix=r1[0]+barcode+'\n+\n'+'I'*12+'\n'
            target.write(ix);wrong.write(r1[0]+'A'*12+'\n+\n'+'I'*12+'\n')
            streams['R1'].writelines(r1);streams['R2'].writelines(r2);streams['I2'].write(ix)
    manifest[number-1].update(index_fastq=str(index),i1_fastq=str(wrong_index))
for stream in streams.values():stream.close()
def write(path,records):
    with path.open('w') as handle:
        writer=csv.DictWriter(handle,list(records[0]));writer.writeheader();writer.writerows(records)
write(fixture/'samples.csv',rows);write(fixture/'fastqs.csv',manifest)
native=[{k:v for k,v in row.items() if k not in ['index','index2','common_linker']} for row in rows]
write(fixture/'native_samples.csv',native)
mock=fixture/'bcl_mock';(mock/'InterOp').mkdir(parents=True,exist_ok=True)
(mock/'InterOp/test.bin').write_bytes(b'routing-test-only')
write(mock/'fastqs.csv',manifest)
write(mock/'samples.csv',[dict(row,Sample_Project='project_'+row['Sample_ID']) for row in rows])
print('Prepared assigned/native Golay FASTQs, empty sample and three-project BCL routing fixture.')
