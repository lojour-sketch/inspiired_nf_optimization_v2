#!/usr/bin/env python3
"""Compare input normalization with the committed upstream metadata helper."""
import importlib.util,json,sys
from pathlib import Path
root=Path(sys.argv[1]).resolve()
spec=importlib.util.spec_from_file_location('inputs',str(root/'bin/bushman_inputs.py'))
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
cases=['N'*12+'CTCCGCTTAAGGGACT',
       'AGCAGGTCCGAAATTCTCGG'+'N'*12+'CTCCGCTTAAGGGACT',
       'ACGTACGTACGTACGTACGT', 'N'*13+'ACGTACGTACGTACGTA']
values=[dict(linkerSequence=x,python_effective=module.original_linker_common(x)) for x in cases]
assert values[0]['python_effective']=='AGTCCCTTAAGCGGAG'
assert values[1]['python_effective']=='AGTCCCTTAAGCGGAG'
for invalid in ['N'*12+'ACGT','ACGT','AAAANAAAANAAAAAAAAAAAAAAAA']:
    try:module.original_linker_common(invalid)
    except ValueError:pass
    else:raise AssertionError('Invalid upstream linker accepted: '+invalid)
Path(sys.argv[2]).write_text(json.dumps(values,indent=2)+'\n')
print('PASS: RUN809 reverse complement and invalid linker layouts')
