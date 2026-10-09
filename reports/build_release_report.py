#!/usr/bin/env python3
"""Build the portable report from release documentation and saved evidence."""
import html,json,re
from pathlib import Path
root=Path(__file__).resolve().parents[1]
def inline(s):
    s=html.escape(s)
    s=re.sub(r'`([^`]+)`',r'<code>\1</code>',s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s)
    return re.sub(r'\[([^]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',s)
def markdown(s):
    out=[];lines=s.splitlines();i=0
    while i<len(lines):
        line=lines[i]
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].startswith('|'):
                row=lines[i]
                if not re.match(r'^\|[\s|:-]+\|$',row):rows.append([inline(x.strip()) for x in row.strip('|').split('|')])
                i+=1
            out.append('<div class="table"><table>'+''.join('<tr>'+''.join('<'+('th' if n==0 else 'td')+'>'+x+'</'+('th' if n==0 else 'td')+'>' for x in row)+'</tr>' for n,row in enumerate(rows))+'</table></div>');continue
        if line.startswith('#'):
            level=len(line)-len(line.lstrip('#'));out.append('<h{0}>{1}</h{0}>'.format(level,inline(line[level:].strip())))
        elif line.startswith('```'):
            code=[];i+=1
            while i<len(lines) and not lines[i].startswith('```'):code.append(lines[i]);i+=1
            out.append('<pre><code>'+html.escape('\n'.join(code))+'</code></pre>')
        elif line.strip():out.append('<p>'+inline(line)+'</p>')
        i+=1
    return '\n'.join(out)
summary=json.loads((root/'reports/validation_summary.json').read_text())
rows=''.join('<tr><td>'+html.escape(item['check'])+'</td><td>'+html.escape(item['status'])+'</td><td>'+html.escape(item['scope'])+'</td></tr>' for item in summary['checks'])
evidence=''.join('<li><a href="evidence/'+p.name+'">'+p.name+'</a></li>' for p in sorted((root/'reports/evidence').glob('*.json')))
links='<p><a href="RUN809_pre_correction/index.html">All seven pre-correction RUN809 plot reports</a> · <a href="../docs/results.md">Results interpretation guide</a></p>'
if (root/'reports/RUN809_corrected_0pct/index.html').exists():links+='<p><a href="RUN809_corrected_0pct/index.html">Corrected 0% validation figures and tables</a></p>'
text='<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Bushman-compatible pipeline: changes and validation</title><style>body{font:16px system-ui;line-height:1.55;margin:40px auto;max-width:1400px;padding:0 24px;color:#1e293b}table{border-collapse:collapse;width:100%}th,td{border:1px solid #cbd5e1;padding:12px;vertical-align:top;text-align:left}th{background:#e2e8f0}.table{overflow-x:auto;margin:24px 0}code{background:#f1f5f9;font-size:.88em;padding:2px}pre{white-space:pre-wrap}.notice{background:#fff3cd;padding:18px}h2{margin-top:38px}a{color:#145ba3}@media print{body{font-size:10pt;padding:0}.table{overflow:visible}}</style></head><body>'
text+='<h1>Bushman-compatible Nextflow pipeline 1.2.0</h1><p>Updated '+html.escape(summary['updated_at_utc'])+'</p><p class="notice">'+html.escape(summary['release_statement'])+'</p>'
text+='<h2>Completed checks and exact validation scope</h2><table><tr><th>Check</th><th>Status</th><th>Scope</th></tr>'+rows+'</table>'+links
text+='<h2>Saved validation evidence</h2><ul>'+evidence+'</ul>'
text+=markdown((root/'docs/objectives_and_structure.md').read_text())
text+=markdown((root/'docs/pipeline_changes.md').read_text())
text+='<h2>Reading the new plots</h2>'+markdown((root/'docs/results.md').read_text())
text+='<h2>Runtime distribution</h2>'+markdown((root/'docs/runtime.md').read_text())+'</body></html>'
(root/'reports/pipeline_changes_report.html').write_text(text)
print('Built comprehensive portable original-Bushman comparison and validation report')
