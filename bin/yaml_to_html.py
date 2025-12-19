#!/usr/bin/env python3
import yaml
import sys
import re
import json

def format_text(text):
    """Format text with markdown-like syntax to HTML with bullet point styling"""
    if not isinstance(text, str):
        return str(text)
    
    # Convert ***bold and italic*** to <strong><em>
    text = re.sub(r'\*\*\*(.+?)\*\*\*', r'<strong><em>\1</em></strong>', text)

    # Convert **bold** to <strong>
    text = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', text)
    
    # Process line by line to handle indentation and bullets
    lines = text.split('\n')
    formatted_lines = []
    
    for line in lines:
        stripped = line.lstrip()
        indent_level = len(line) - len(stripped)
        # Empty lines create a line break if not stripped: 
        if line.strip() == '':
            formatted_lines.append('<br>') 
            continue
        
        # Check if line starts with bullet
        if stripped.startswith('- '):
            bullet_text = stripped[2:]  # Remove "- "
            
            # Different bullet styles based on indentation
            if indent_level == 3:
                # Main level - filled circle
                formatted_lines.append(f'<span style="display:block; margin-left:{indent_level}px;">➛ {bullet_text}</span>')
            elif indent_level == 4:
                # First sub-level - empty circle
                formatted_lines.append(f'<span style="display:block; margin-left:{indent_level*2}px;">• {bullet_text}</span>')
            elif indent_level == 5:
                # Second sub-level - square
                formatted_lines.append(f'<span style="display:block; margin-left:{indent_level*3}px;">▪ {bullet_text}</span>')
            else:
                # Third+ sub-level - dash
                formatted_lines.append(f'<span style="display:block; margin-left:{indent_level*4}px;">✧ {bullet_text}</span>')
        else:
            # Regular text with indentation
            if indent_level > 0:
                formatted_lines.append(f'<span style="display:block; margin-left:{indent_level*4}px;">{stripped}</span>')
            else:
                formatted_lines.append(stripped)
    
    return ''.join(formatted_lines)

def yaml_to_html(yaml_file, output_file):
    """Convert meta.yml to readable HTML"""
    
    try:
        with open(yaml_file, 'r') as f:
            data = yaml.safe_load(f)
    except Exception as e:
        print(f"❌ Error parsing YAML: {e}")
        sys.exit(1)
    
    # Debug: Print structure
    print("🔍 Debugging YAML structure:")
    print(f"   Keys found: {list(data.keys())}")
    if 'output' in data:
        print(f"   Output type: {type(data['output'])}")
        print(f"   Output content preview: {json.dumps(data['output'], indent=2, default=str)[:500]}")
    
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Module Documentation</title>
    <style>
        * { box-sizing: border-box; }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            max-width: 1400px;
            margin: 0 auto;
            padding: 40px 20px;
            background-color: #f8f9fa;
            color: #212529;
            line-height: 1.7;
        }
        
        .container {
            background: white;
            padding: 40px;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        
        h1 {
            color: #24B064;
            font-size: 2.5em;
            margin-bottom: 10px;
            border-bottom: 4px solid #24B064;
            padding-bottom: 15px;
        }
        
        h2 {
            color: #2c3e50;
            font-size: 1.8em;
            margin-top: 40px;
            margin-bottom: 20px;
            border-bottom: 2px solid #e9ecef;
            padding-bottom: 10px;
        }
        
        h3 {
            color: #495057;
            font-size: 1.3em;
            margin-top: 25px;
            margin-bottom: 15px;
            margin-left: 0;
        }
        
        .section {
            margin: 30px 0;
        }
        
        .description {
            background-color: #f8f9fa;
            padding: 25px;
            border-left: 5px solid #24B064;
            margin: 25px 0;
            font-size: 0.95em;
            line-height: 1.8;
        }
        
        .tool-info {
            background-color: #f8f9fa;
            padding: 15px 20px;
            margin: 15px 0 15px 20px;
            border-radius: 5px;
            border-left: 3px solid #6c757d;
        }
        
        .tool-info p {
            margin: 8px 0;
        }
        
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.92em;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
        }
        
        th {
            background-color: #24B064;
            color: white;
            font-weight: 600;
            text-align: left;
            padding: 15px 12px;
            text-transform: uppercase;
            font-size: 0.85em;
            letter-spacing: 0.5px;
        }
        
        td {
            padding: 15px 12px;
            border-bottom: 1px solid #e9ecef;
            vertical-align: top;
        }
        
        tr:hover {
            background-color: #f8f9fa;
        }
        
        tr:last-child td {
            border-bottom: none;
        }
        
        code {
            background-color: #e9ecef;
            padding: 3px 8px;
            border-radius: 4px;
            font-family: 'Consolas', 'Monaco', 'Courier New', monospace;
            font-size: 0.9em;
            color: #e83e8c;
            font-weight: 500;
        }
        
        .pattern {
            color: #6c757d;
            font-style: italic;
            font-size: 0.9em;
        }
        
        ul {
            list-style-type: disc;
            padding-left: 25px;
            margin: 15px 0;
        }
        
        li {
            padding: 5px 0;
            line-height: 1.6;
        }
        
        a {
            color: #007bff;
            text-decoration: none;
            border-bottom: 1px solid transparent;
            transition: border-color 0.2s;
        }
        
        a:hover {
            border-bottom-color: #007bff;
        }
        
        .keywords {
            background-color: #e7f3ff;
            padding: 12px 20px;
            border-radius: 5px;
            margin: 20px 0;
        }
        
        .keywords strong {
            color: #0056b3;
        }
        
        .note {
            background-color: #fff3cd;
            border-left: 4px solid #ffc107;
            padding: 15px 20px;
            margin: 20px 0;
            border-radius: 4px;
        }
        
        .debug {
            background-color: #f0f0f0;
            padding: 10px;
            margin: 10px 0;
            border: 1px solid #ccc;
            font-family: monospace;
            font-size: 0.85em;
        }
    </style>
</head>
<body>
    <div class="container">
"""
    
    # Title
    html_content += f"        <h1>{data.get('name', 'Module Documentation')}</h1>\n"
    
    # Description
    if 'description' in data:
        formatted_desc = format_text(data['description'])
        html_content += f'        <div class="description">{formatted_desc}</div>\n'
    
    # Keywords
    if 'keywords' in data:
        keywords_str = ', '.join(data['keywords'])
        html_content += f'        <div class="keywords"><strong>Keywords:</strong> {keywords_str}</div>\n'
    
    # Tools
    if 'tools' in data:
        html_content += '        <div class="section">\n            <h2>Tools</h2>\n'
        for tool in data['tools']:
            for tool_name, tool_info in tool.items():
                html_content += f'            <h3>{tool_name}</h3>\n'
                html_content += '            <div class="tool-info">\n'
                html_content += f'                <p>{tool_info.get("description", "")}</p>\n'
                if 'homepage' in tool_info:
                    html_content += f'                <p><strong>Homepage:</strong> <a href="{tool_info["homepage"]}" target="_blank">{tool_info["homepage"]}</a></p>\n'
                if 'documentation' in tool_info:
                    html_content += f'                <p><strong>Documentation:</strong> <a href="{tool_info["documentation"]}" target="_blank">{tool_info["documentation"]}</a></p>\n'
                html_content += '            </div>\n'
        html_content += '        </div>\n'
    
    # Input
    if 'input' in data:
        html_content += '        <div class="section">\n            <h2>Input Parameters</h2>\n'
        html_content += '            <table>\n'
        html_content += '                <tr><th>Name</th><th>Type</th><th>Description</th><th>Pattern</th></tr>\n'
        
        for inp in data['input']:
            if isinstance(inp, list):
                for item in inp:
                    if isinstance(item, dict):
                        for name, details in item.items():
                            if isinstance(details, dict):
                                pattern = details.get('pattern', '')
                                inp_type = details.get('type', '')
                                description = format_text(details.get('description', ''))
                                html_content += f'                <tr><td><code>{name}</code></td><td>{inp_type}</td><td>{description}</td><td class="pattern">{pattern}</td></tr>\n'
            elif isinstance(inp, dict):
                for name, details in inp.items():
                    if isinstance(details, dict):
                        pattern = details.get('pattern', '')
                        inp_type = details.get('type', '')
                        description = format_text(details.get('description', ''))
                        html_content += f'                <tr><td><code>{name}</code></td><td>{inp_type}</td><td>{description}</td><td class="pattern">{pattern}</td></tr>\n'
        
        html_content += '            </table>\n        </div>\n'
    
    # Output - ONE ROW PER CHANNEL WITH ALL ELEMENTS
    if 'output' in data and data['output']:
        html_content += '        <div class="section">\n            <h2>Output Channels</h2>\n'
        
        html_content += '            <table>\n'
        html_content += '                <tr><th>Channel</th><th>Element</th><th>Type</th><th>Description</th><th>Pattern</th></tr>\n'
        
        row_count = 0
        
        # Process each output channel
        for out in data['output']:
            if isinstance(out, dict):
                for channel_name, channel_details in out.items():
                    
                    # If channel_details is a list (e.g., allsites: [list of file dicts])
                    if isinstance(channel_details, list):
                        first_element = True
                        
                        for item in channel_details:
                            if isinstance(item, dict):
                                # Each dict in the list is a file
                                for file_name, file_details in item.items():
                                    channel_cell = f'<code>{channel_name}</code>' if first_element else ''
                                    
                                    if isinstance(file_details, dict):
                                        pattern = file_details.get('pattern', '')
                                        file_type = file_details.get('type', 'file')
                                        description = format_text(file_details.get('description', ''))
                                        
                                        html_content += f'                <tr><td>{channel_cell}</td><td><code>{file_name}</code></td><td>{file_type}</td><td>{description}</td><td class="pattern">{pattern}</td></tr>\n'
                                    else:
                                        # file_details is just a pattern string
                                        html_content += f'                <tr><td>{channel_cell}</td><td><code>{file_name}</code></td><td>file</td><td></td><td class="pattern">{file_details}</td></tr>\n'
                                    
                                    first_element = False
                                    row_count += 1
                    
                    # If channel_details is a dict (single file output)
                    elif isinstance(channel_details, dict):
                        pattern = channel_details.get('pattern', '')
                        file_type = channel_details.get('type', 'file')
                        description = format_text(channel_details.get('description', ''))
                        
                        html_content += f'                <tr><td><code>{channel_name}</code></td><td>-</td><td>{file_type}</td><td>{description}</td><td class="pattern">{pattern}</td></tr>\n'
                        row_count += 1
        
        if row_count == 0:
            html_content += '                <tr><td colspan="5" style="text-align:center; color:#999;">No output channels could be parsed. Check YAML structure.</td></tr>\n'
        
        html_content += '            </table>\n        </div>\n'
    else:
        html_content += '        <div class="note"><strong>⚠️ Note:</strong> No output section found in meta.yml</div>\n'

    # Authors
    if 'authors' in data:
        html_content += '        <div class="section">\n            <h2>Authors</h2>\n            <ul>\n'
        authors = data['authors']
        if isinstance(authors, list):
            for author in authors:
                html_content += f'                <li>{str(author)}</li>\n'
        else:
            html_content += f'                <li>{str(authors)}</li>\n'
        html_content += '            </ul>\n        </div>\n'
    
    # Maintainers
    if 'maintainers' in data:
        html_content += '        <div class="section">\n            <h2>Maintainers</h2>\n            <ul>\n'
        maintainers = data['maintainers']
        if isinstance(maintainers, list):
            for maintainer in maintainers:
                html_content += f'                <li>{str(maintainer)}</li>\n'
        else:
            html_content += f'                <li>{str(maintainers)}</li>\n'
        html_content += '            </ul>\n        </div>\n'
    
    html_content += """    </div>
</body>
</html>"""
    
    with open(output_file, 'w') as f:
        f.write(html_content)
    
    print(f"\n✅ Successfully created {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python yaml_to_html_final.py <meta.yml> <output.html>")
        print("Example: python yaml_to_html_final.py modules/local/alignment/meta.yml docs/alignment.html")
        sys.exit(1)
    
    yaml_to_html(sys.argv[1], sys.argv[2])