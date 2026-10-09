#!/usr/bin/env python3
"""Run the original Golay decoder with Python's eager-map compatibility fix."""
import sys
from pathlib import Path
import types
sys.path.insert(0,'/opt/bushman_python' if Path('/opt/bushman_python').is_dir() else
    str(Path(__file__).resolve().parent.parent/'vendor/bushman_python'))
code=Path(sys.argv[1]).read_text().replace('numpy.array(map(int,bitstring))','numpy.array(list(map(int,bitstring)))')
module=types.ModuleType('original_golay')
exec(compile(code,sys.argv[1],'exec'),module.__dict__)
for path in sys.argv[2:]:
    output=Path(path).with_name(Path(path).name.replace('trimmedI1','correctedI1'))
    with open(path) as inp,output.open('w') as out:
        name=None;parts=[]
        def write():
            if name is not None:
                corrected=module.decode(''.join(parts))[0]
                if corrected is not None:out.write('>'+name+'\n'+corrected+'\n')
        for line in inp:
            if line.startswith('>'):write();name=line[1:].split()[0];parts=[]
            else:parts.append(line.strip())
        write()
